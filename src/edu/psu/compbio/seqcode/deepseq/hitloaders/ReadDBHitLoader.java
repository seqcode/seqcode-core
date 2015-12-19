package edu.psu.compbio.seqcode.deepseq.hitloaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import edu.psu.compbio.seqcode.deepseq.HitPair;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.projects.readdb.*;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * ReadDBHitLoader: load read alignment five primes from a collection of ReadDB alignments (SeqLocators)
 * @author mahony
 *
 */
public class ReadDBHitLoader extends HitLoader{
	
	private Genome gen=null;
	private Client client=null; //ReadDB client
	private SeqDataLoader loader=null; 
	private List<SeqLocator> exptLocs;
	private List<String> exptNames =new ArrayList<String>();
	private List<SeqAlignment> aligns = new ArrayList<SeqAlignment>();
	private Collection<String> alignIDs= new ArrayList<String>();
	private HashMap<SeqAlignment, Set<Integer>> availSingleChroms = new HashMap<SeqAlignment, Set<Integer>>();
	private HashMap<SeqAlignment, Set<Integer>> availSingleType2Chroms = new HashMap<SeqAlignment, Set<Integer>>();
	private HashMap<SeqAlignment, Set<Integer>> availPairedChroms = new HashMap<SeqAlignment, Set<Integer>>();
	private boolean hasPairedAligns=false;
	private final int MAXRDBLOAD = 1000000; //Maximum number of hits to load from ReadDB in one chunk
	
	/**
	 * Constructor: initialize the experiments.
	 * @param g Genome
	 * @param locs ChipSeqLocator
	 * @param loadType1Reads boolean (load standard single-end/paired-end reads)
	 * @param loadType2Reads boolean (load type 2 reads if they exist)
	 * @param loadPairs boolean
	 */
	public ReadDBHitLoader(SeqDataLoader loader, Genome g, List<SeqLocator> locs, boolean loadT1Reads, boolean loadT2Reads, boolean loadPairs){
		super(loadT1Reads, loadT2Reads, true, loadPairs);
		gen=g;
		exptLocs = locs;
		
		if(gen==null){
			System.err.println("ReadDBHitLoader: null genome provided. ReadDB requires a defined genome. ");
		}
		if(!loadType1 && !loadType2 && !loadPairs){
			System.err.println("ReadDBHitLoader: Error: you didn't request to load any read types. ");
		}
        if (exptLocs.size() == 0) {
            System.err.println("Created a ReadDBHitLoader with no SeqLocators");
        }

		try {
			//Start a new ReadDB client
			if(client==null)
				client = new Client();

			//Process the SeqLocators
			for(SeqLocator locator : exptLocs){
				String exptName = locator.getExptName(); exptNames.add(exptName);
				aligns.addAll(loader.loadAlignments(locator, gen));
			}
	        for(SeqAlignment alignment : aligns) {
	        	if(client.exists(Integer.toString(alignment.getDBID()))){
		            alignIDs.add(Integer.toString(alignment.getDBID()));
		            this.sourceName=this.sourceName.equals("") ? 
		            		alignment.getExpt().getName()+";"+alignment.getExpt().getReplicate()+";"+alignment.getName() :
		            		this.sourceName+":"+alignment.getExpt().getName()+";"+alignment.getExpt().getReplicate()+";"+alignment.getName();
		            hasPairedAligns = hasPairedAligns || (alignment.getAlignType().getName().equals("PAIRED") || alignment.getAlignType().getName().equals("MIXED"));
		            
		            //Find the available chromosomes for each alignment
		            availSingleChroms.put(alignment, new HashSet<Integer>());
		            if(alignment.getNumHits()>0){
		            	availSingleChroms.get(alignment).addAll(client.getChroms(Integer.toString(alignment.getDBID()), false,false, null));
		            }
		        	availSingleType2Chroms.put(alignment, new HashSet<Integer>());
		            if(alignment.getNumType2Hits()>0){
		            	availSingleType2Chroms.get(alignment).addAll(client.getChroms(Integer.toString(alignment.getDBID()), true,false, null));
		            }
		            availPairedChroms.put(alignment, new HashSet<Integer>());
		            if(alignment.getNumPairs()>0){
		            	availPairedChroms.get(alignment).addAll(client.getChroms(Integer.toString(alignment.getDBID()), false,true, null));
		            }
	        	}else{
	        		System.err.println("ReadDBHitLoader: Error: "+alignment.getExpt().getName()+";"+alignment.getExpt().getReplicate()+";"+alignment.getName()+"\tRDBID:"+alignment.getDBID()+" does not exist in ReadDB.");
	        		System.exit(1);
	        	}
		    }
	        
            
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClientException e) {
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
	public void sourceAllHits(){
		this.initialize();
		try {
			//Start a new ReadDB client
			if(client==null)
				client = new Client();

			//Iterate over each chromosome
			for (String chrom: gen.getChromList()){
				// load  data for this chromosome.
				int length = gen.getChromLength(chrom);
				Region wholeChrom = new Region(gen, chrom, 1, length);
				int count = 0;
				for(SeqAlignment alignment : aligns) { 
					if(loadType1)
						if(availSingleChroms.get(alignment).contains(gen.getChromID(wholeChrom.getChrom()))){
							count += client.getCount(Integer.toString(alignment.getDBID()),
		                				gen.getChromID(wholeChrom.getChrom()),
		                                false,
		                				false,
		                                wholeChrom.getStart(),
		                                wholeChrom.getEnd(),
		                                null,
		                                null,
		                                null);
						}
					if(loadType2){
						if(availSingleType2Chroms.get(alignment).contains(gen.getChromID(wholeChrom.getChrom()))){
							count += client.getCount(Integer.toString(alignment.getDBID()),
		                				gen.getChromID(wholeChrom.getChrom()),
		                                true,
		                				false,
		                                wholeChrom.getStart(),
		                                wholeChrom.getEnd(),
		                                null,
		                                null,
		                                null);
						}
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
					Pair<ArrayList<Integer>,ArrayList<Float>> hits;
					if(loadType1){
						hits = loadStrandedBaseCounts(chunk, '+', false);
						addHits(chrom, '+', hits.car(), hits.cdr());
						hits = loadStrandedBaseCounts(chunk, '-', false);
						addHits(chrom, '-', hits.car(), hits.cdr());
					}
					if(loadType2){
						hits = loadStrandedBaseCounts(chunk, '+', true);
						addHits(chrom, '+', hits.car(), hits.cdr());
						hits = loadStrandedBaseCounts(chunk, '-', true);
						addHits(chrom, '-', hits.car(), hits.cdr());
					}
					
					if(hasPairedAligns && loadPairs){
						ArrayList<HitPair> pairs = loadStrandedPairs(chunk, '+');
						addPairs(chrom, '+', pairs);
						pairs = loadStrandedPairs(chunk, '-');
						addPairs(chrom, '-', pairs);
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here: ClientException could be thrown because chromosome doesn't contain any hist. 
		}
	}
	
	
    /**
     *  load pairs of read hit 5' coordinates (sorted) and counts
     * 
     */
    private Pair<ArrayList<Integer>,ArrayList<Float>> loadStrandedBaseCounts(Region r, char strand, boolean loadRead2){
    	
        TreeMap<Integer,Float> allHits = null;
        ArrayList<Integer> coords = new ArrayList<Integer>();
        ArrayList<Float> counts = new ArrayList<Float>();
        try {
        	//Start a new ReadDB client
    		if(client==null)
    			client = new Client();

    		allHits = client.getWeightHistogram(alignIDs,
                                                r.getGenome().getChromID(r.getChrom()),
                                                loadRead2,
                                                false,
                                                0,
                                                1,
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                strand == '+');
            
            if (allHits == null) {
                if (alignIDs.size() != 0) {
                    //throw new NullPointerException("ReadDBHitLoader: client.getWeightHistogram returned null");
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
     *  load paired reads in region
     * 
     */
    private ArrayList<HitPair> loadStrandedPairs(Region r, char strand){
    	ArrayList<HitPair> pairs = new ArrayList<HitPair>();
    	try{
    		//Start a new ReadDB client
    		if(client==null)
    			client = new Client();
    		
	    	for (String alignid : alignIDs) {
	            List<PairedHit> ph = client.getPairedHits(alignid,
	                                                     r.getGenome().getChromID(r.getChrom()),
	                                                     true,
	                                                     r.getStart(),
	                                                     r.getEnd(),
	                                                     null,
	                                                     strand=='+');
	            for (PairedHit h : ph) {
	            	if(h.pairCode==1)
	            		pairs.add(new HitPair(h.leftPos, r.getGenome().getChromName(h.rightChrom), h.rightPos, h.rightStrand?0:1, h.weight));
	            }
	    	}
    	} catch (IOException e) {
            e.printStackTrace();
        } catch (ClientException e) {
            //Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
        }
    	return pairs;
    }
	
    /**
     * Close the client
     */
	public void cleanup(){
		if(client!=null)
			client.close();
		if(loader!=null)
			loader.close();
	}
}
