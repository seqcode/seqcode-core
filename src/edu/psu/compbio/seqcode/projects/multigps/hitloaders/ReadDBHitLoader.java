package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.projects.readdb.*;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * ReadDBHitLoader: load read alignment five primes from a collection of ReadDB alignments (ChipSeqLocators)
 * @author mahony
 *
 */
public class ReadDBHitLoader extends HitLoader{
	
	private Genome gen=null;
	private Client client=null; //ReadDB client
	private List<SeqLocator> exptLocs;
	private List<String> exptNames =new ArrayList<String>();
	private List<SeqAlignment> aligns = new ArrayList<SeqAlignment>();
	private Collection<String> alignIDs= new ArrayList<String>();
	private final int MAXRDBLOAD = 1000000; //Maximum number of hits to load from ReadDB in one chunk
	
	/**
	 * Constructor: initialize the experiments
	 * @param g Genome
	 * @param locs ChipSeqLocator
	 */
	public ReadDBHitLoader(Genome g, List<SeqLocator> locs){
		super();
		gen=g;
		exptLocs = locs;
		if(gen==null){
			System.err.println("ReadDBHitLoader: null genome provided. ReadDB requires a defined genome. ");
		}
        if (exptLocs.size() == 0) {
            System.err.println("Created a ReadDBHitLoader with no ChipSeqLocators");
        }

		try {
			
			//Initialize ChipSeqLoaders
            SeqDataLoader loader = new SeqDataLoader(false); 
			for(SeqLocator locator : exptLocs){
				String exptName = locator.getExptName(); exptNames.add(exptName);
				if (locator.getAlignName() == null) {
		            if(locator.getReplicates().isEmpty()) { //No alignment name, no replicate names
		            	Collection<SeqExpt> expts = loader.loadExperiments(locator.getExptName());
		        		for(SeqExpt expt : expts) { 
		                	Collection<SeqAlignment> aligns;
							aligns = loader.loadAllAlignments(expt);
							for (SeqAlignment currentAlign : aligns) {
		            			if (currentAlign.getGenome().equals(g)) { 
		            				aligns.add(currentAlign);
		    						break;
		    					}
		            		}
		    			}
		            } else { //No alignment name, given replicate names
		                for(String repName : locator.getReplicates()) { 
		                    SeqExpt expt = loader.loadExperiment(locator.getExptName(), repName);
		                    SeqAlignment alignment = 
		                        loader.loadAlignment(expt, locator.getAlignName(), g);
		                    if(alignment != null) { 
		                        aligns.add(alignment);
		                        break;
		                    }
		                }
		            }
		        } else {
		        	if(locator.getReplicates().isEmpty()) {//Given alignment name, no replicate names
		        		Collection<SeqExpt> expts = loader.loadExperiments(locator.getExptName());
                        System.err.println("Have name but no replicates.  Got " + expts.size() + " experiments for " + locator.getExptName());
		        		for(SeqExpt expt : expts) { 
		                	Collection<SeqAlignment> alignments;
							alignments = loader.loadAllAlignments(expt);
							for (SeqAlignment currentAlign : alignments) {
                                System.err.println("  " + currentAlign);
		            			if (currentAlign.getGenome().equals(g) && currentAlign.getName().equals(locator.getAlignName())) { 
		            				aligns.add(currentAlign);
		    						break;
		    					}
		            		}
		    			}
		            }else{
		            	for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
		            		SeqAlignment a = loader.loadAlignment(loader.loadExperiment(locator.getExptName(),replicate),locator.getAlignName(),g);
							if(a!=null)
								aligns.add(a);
		        		}
		            }
		        }
			}
	        for(SeqAlignment alignment : aligns) {
	            alignIDs.add(Integer.toString(alignment.getDBID()));
	        }
	        
            if (exptLocs.size() != 0 && aligns.size() == 0) {
                System.err.println("Locators were " + exptLocs + " but didn't get any alignments");
            }
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * Load the five primes from ReadDB
	 */
	public void sourceReads(){
		this.initialize();
		try {
			//Start a new ReadDB client
			if(client==null)
				client = new Client();
			
			//Find the available chromosomes for each alignment
			HashMap<SeqAlignment, Set<Integer>> availChroms = new HashMap<SeqAlignment, Set<Integer>>();
			for(SeqAlignment alignment : aligns) {
				availChroms.put(alignment, client.getChroms(Integer.toString(alignment.getDBID()), false, null));
			}
			
			//Iterate over each chromosome
			for (String chrom: gen.getChromList()){
				// load  data for this chromosome.
				int length = gen.getChromLength(chrom);
				Region wholeChrom = new Region(gen, chrom, 1, length);
				int count = 0;
				for(SeqAlignment alignment : aligns) { 
					if(availChroms.get(alignment).contains(gen.getChromID(wholeChrom.getChrom()))){
		                count += client.getCount(Integer.toString(alignment.getDBID()),
		                				gen.getChromID(wholeChrom.getChrom()),
		                                false,
		                                wholeChrom.getStart(),
		                                wholeChrom.getEnd(),
		                                null,
		                                null,
		                                null);
					}
				}
				ArrayList<Region> chunks = new ArrayList<Region>();
				// if there are too many reads in a chrom, read smaller chunks
				if (count>MAXRDBLOAD){
					int chunkNum = count/MAXRDBLOAD*2+1;
					int chunkLength = length/chunkNum;
					int start = 0;
					while (start<=length){
						int end = Math.min(length, start+chunkLength-1);
						Region r = new Region(gen, chrom, start, end);
						start = end+1;
						chunks.add(r);
					}
				}else
					chunks.add(wholeChrom);

				for (Region chunk: chunks){
					Pair<ArrayList<Integer>,ArrayList<Float>> hits = loadStrandedBaseCounts(chunk, '+');
					addHits(chrom, '+', hits.car(), hits.cdr());
					hits = loadStrandedBaseCounts(chunk, '-');
					addHits(chrom, '-', hits.car(), hits.cdr());
				}
			}
			
			client.close();
			client=null;
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClientException e) {
			e.printStackTrace();
		}
	}
	
    /**
     *  load paired read hit 5' coordinates (sorted) and counts
     * 
     */
    private Pair<ArrayList<Integer>,ArrayList<Float>> loadStrandedBaseCounts(Region r, char strand){
        
        TreeMap<Integer,Float> allHits = null;
        ArrayList<Integer> coords = new ArrayList<Integer>();
        ArrayList<Float> counts = new ArrayList<Float>();
        try {
            allHits = client.getWeightHistogram(alignIDs,
                                                r.getGenome().getChromID(r.getChrom()),
                                                false,
                                                false,
                                                1,
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                strand == '+');
            if (allHits == null) {
                if (alignIDs.size() != 0) {
                    throw new NullPointerException("ReadDBHitLoader: client.getWeightHistogram returned null");
                }
            } else {
                coords.addAll(allHits.keySet());
                counts.addAll(allHits.values());
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClientException e) {
            //Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
        }
        return new Pair<ArrayList<Integer>,ArrayList<Float>>(coords, counts);
    }
	
    /**
     * Close the client
     */
	public void cleanup(){
		if(client!=null){
			client.close();
		}
	}
}
