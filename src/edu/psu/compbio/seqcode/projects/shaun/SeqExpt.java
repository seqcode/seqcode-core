package edu.psu.compbio.seqcode.projects.shaun;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import edu.psu.compbio.seqcode.gse.datasets.chippet.RunningOverlapSum;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqHit;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLoader;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.ChipSeqExpander;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqExpt {
    public static final int defaultReadLength = 26, defaultReadExtension = 174;
	private double readLength = defaultReadLength;
	private double readExtension = defaultReadExtension;
	private String exptName;
	private Genome currentGen=null;
	private double hitCount = 0;
	private double totalSeq = 0;
	private double expect=0;
	private ChipSeqExpander expander=null;
	private ChipSeqLocator loc = null;
	private ChipSeqLoader loader;
	private LinkedList<ChipSeqAlignment> alignments;

	public SeqExpt(Genome g, String exptName){
		this(g, exptName, "");
	}
	 
	public SeqExpt(Genome g,  String exptName, String replicate){
		currentGen=g;
		this.exptName=exptName;
		
		ChipSeqLoader chipSeqLoader = null;
    	try {
    		chipSeqLoader = new ChipSeqLoader();
    	} 
    	catch (SQLException sqlex) {
    		sqlex.printStackTrace();
    	} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
        if (chipSeqLoader != null) { 
        	try{        		
        		if (replicate.equals("")) {
        			List<ChipSeqLocator> locs = new Vector<ChipSeqLocator>();
        			
        			List<ChipSeqExpt> expts = new Vector<ChipSeqExpt>();
        			expts.addAll(chipSeqLoader.loadExperiments(exptName));
					
        			for (ChipSeqExpt expt : expts) {
        				Collection<ChipSeqAlignment> aligns;
						aligns = chipSeqLoader.loadAllAlignments(expt);
						
                		for (ChipSeqAlignment currentAlign : aligns) {
                			if (currentAlign.getGenome().equals(g)) { 
                				ChipSeqLocator currentLoc = new ChipSeqLocator(expt.getName(), 
                                        expt.getReplicate(), currentAlign.getName());
                				locs.add(currentLoc);
        						break;
        					}
                		}
        			}

        			List<ChipSeqLocator> collapsedLocs = new Vector<ChipSeqLocator>(this.collapseLocatorsByName(locs));
        			if (collapsedLocs.size() != 1) {
        				System.err.println(collapsedLocs.size() + " collapsed locators");
        				System.exit(0);
        			}
        			loc = collapsedLocs.get(0);
        		}
        		else {
        			ChipSeqExpt expt;
					expt = chipSeqLoader.loadExperiment(exptName, replicate);
					
        			ChipSeqAlignment align = null;
            		Collection<ChipSeqAlignment> aligns;
					aligns = chipSeqLoader.loadAllAlignments(expt);
					
            		for (ChipSeqAlignment currentAlign : aligns) {
            			if (currentAlign.getGenome().equals(g)) { 
    						align = currentAlign;
    						break;
    					}
            		}
            		
            		loc = new ChipSeqLocator(expt.getName(), expt.getReplicate(), align.getName());
        		}
        		
        		expander = new ChipSeqExpander(loc);
        	    loader = new ChipSeqLoader();
    	        alignments = new LinkedList<ChipSeqAlignment>();
    	        
    	        if(loc.getReplicates().isEmpty()) { 
    	        	Collection<ChipSeqExpt> expts = loader.loadExperiments(loc.getExptName());
    	        	for(ChipSeqExpt expt : expts) { 
    	        		ChipSeqAlignment alignment = 
    	        			loader.loadAlignment(expt, loc.getAlignName(), g);
    	        		if(alignment != null) { 
    	        			alignments.add(alignment);
    	        		}
    	        	}
    	        } else {
    	            for(String repName : loc.getReplicates()) { 
    	                ChipSeqExpt expt = loader.loadExperiment(loc.getExptName(), repName);
    	                ChipSeqAlignment alignment = 
    	                    loader.loadAlignment(expt, loc.getAlignName(), g);
    	                if(alignment != null) { 
    	                    alignments.add(alignment);
    	                }
    	            }
    	        }
    	        countHits();
	        } catch (SQLException e) {
				e.printStackTrace();
			} catch (NotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
        }        
	}
    public SeqExpt(Genome g, ChipSeqLocator locator) throws NotFoundException, SQLException {
		currentGen=g;
        exptName = locator.getExptName();
        try {
			loader = new ChipSeqLoader();
			
	        alignments = new LinkedList<ChipSeqAlignment>();
	        if (locator.getAlignName() == null) {
	            if(locator.getReplicates().isEmpty()) { //No alignment name, no replicate names
	            	Collection<ChipSeqExpt> expts = loader.loadExperiments(locator.getExptName());
	        		List<ChipSeqLocator> locs = new Vector<ChipSeqLocator>();
	        		for(ChipSeqExpt expt : expts) { 
	                	Collection<ChipSeqAlignment> aligns;
						aligns = loader.loadAllAlignments(expt);
						for (ChipSeqAlignment currentAlign : aligns) {
	            			if (currentAlign.getGenome().equals(g)) { 
	            				ChipSeqLocator currentLoc = new ChipSeqLocator(expt.getName(), 
	                                    expt.getReplicate(), currentAlign.getName());
	            				locs.add(currentLoc);
	            				alignments.add(currentAlign);
	    						break;
	    					}
	            		}
	    			}
	    			List<ChipSeqLocator> collapsedLocs = new Vector<ChipSeqLocator>(this.collapseLocatorsByName(locs));
	    			if (collapsedLocs.size() != 1) {
	    				System.err.println(collapsedLocs.size() + " collapsed locators");
	    				System.exit(0);
	    			}
	    			locator = collapsedLocs.get(0);
	            } else { //No alignment name, given replicate names
	                for(String repName : locator.getReplicates()) { 
	                    ChipSeqExpt expt = loader.loadExperiment(locator.getExptName(), repName);
	                    ChipSeqAlignment alignment = 
	                        loader.loadAlignment(expt, locator.getAlignName(), g);
	                    if(alignment != null) { 
	                        locator = new ChipSeqLocator(locator.getExptName(),
	                                                     locator.getReplicates(),
	                                                     alignment.getName());
	                        alignments.add(alignment);
	                        break;
	                    }
	                }
	            }
	        } else {
	        	if(locator.getReplicates().isEmpty()) {//Given alignment name, no replicate names
	        		Collection<ChipSeqExpt> expts = loader.loadExperiments(locator.getExptName());
	        		List<ChipSeqLocator> locs = new Vector<ChipSeqLocator>();
	        		for(ChipSeqExpt expt : expts) { 
	                	Collection<ChipSeqAlignment> aligns;
						aligns = loader.loadAllAlignments(expt);
						for (ChipSeqAlignment currentAlign : aligns) {
	            			if (currentAlign.getGenome().equals(g) && currentAlign.getName().equals(locator.getAlignName())) { 
	            				ChipSeqLocator currentLoc = new ChipSeqLocator(expt.getName(), 
	                                    expt.getReplicate(), currentAlign.getName());
	            				locs.add(currentLoc);
	            				alignments.add(currentAlign);
	    						break;
	    					}
	            		}
	    			}
	    			List<ChipSeqLocator> collapsedLocs = new Vector<ChipSeqLocator>(this.collapseLocatorsByName(locs));
	    			if (collapsedLocs.size() != 1) {
	    				System.err.println(collapsedLocs.size() + " collapsed locators");
	    				System.exit(0);
	    			}
	    			locator = collapsedLocs.get(0);
	
	            }else{
	            	for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
	        			alignments.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(),
                                                                                  replicate), 
                                                            locator.getAlignName(),
                                                            g));
	        		}
	            }
	        }


	        expander = new ChipSeqExpander(locator);
	        loc = locator;
	        countHits();
        } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
    }
	
	public double getHitCount(){return(hitCount);}
	public double getTotalSeq(){return(totalSeq);}
	public void setReadLength(double rl){readLength=rl;}
	public void setReadExtension(double re){readExtension=re;}
	public double getReadLength(){return(readLength);}
	public double getExtendedReadLength(){return(readLength+readExtension);}
	public double getExtension(){return readExtension;}
	public Genome getGenome(){return currentGen;}
	public String getExptName(){return exptName;}
	public ChipSeqExpander getExpander(){return expander;}
	public ChipSeqLocator getLocator(){return loc;}
	
	
	public double countHits(){
		hitCount=0;
		totalSeq=0;
		Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			
			hitCount += expander.getHitCount(currentChrom);
			totalSeq += (double)currentChrom.getWidth();
		}		
		
		if(hitCount >0){
			expect = ((readLength+readExtension)*(hitCount/totalSeq));
		}
		
		return hitCount;
	}
	public double countHits(Region a) {
		double hits=0;
		try {
			for(ChipSeqAlignment alignment : alignments) { 
				double currHits = (double)loader.countByRegion(alignment, a);
				hits+=currHits;
			}
			return hits;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0;
		}
	}
	
	public LinkedList<StrandedRegion> loadHits(Region a) {
		try {
			LinkedList<StrandedRegion> total = new LinkedList<StrandedRegion>();
			for(ChipSeqAlignment alignment : alignments) { 
				Collection<ChipSeqHit> hits = loader.loadByRegion(alignment, a);
				for(ChipSeqHit curr : hits){
					total.add(hit2region(0, curr));
				}
			}
			return total;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return new LinkedList<StrandedRegion>();
		}
	}
	
	public LinkedList<StrandedRegion> loadExtendedHits(Region a) {
		try {
			LinkedList<StrandedRegion> total = new LinkedList<StrandedRegion>();
			for(ChipSeqAlignment alignment : alignments) { 
				Collection<ChipSeqHit> hits = loader.loadByRegion(alignment, a);
				for(ChipSeqHit curr : hits){
					total.add(hit2region((int)readExtension, curr));
				}
			}
			return total;
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return new LinkedList<StrandedRegion>();
		}
	}
	public int [] hitCountLandscape(Region a) {
		int [] hits = new int [a.getWidth()+1];
		for(int i=0; i<=a.getWidth(); i++){
			hits[i]=0;
		}
		
		try {
			for(ChipSeqAlignment alignment : alignments) { 
				Iterator<ChipSeqHit> h = loader.loadByRegion(alignment, a).iterator();
				while(h.hasNext()){
					ChipSeqHit csh = h.next();
					int pos = csh.getStart()-a.getStart();
					if(pos<0){pos=0;}
					if(pos>a.getWidth())pos = a.getWidth();
					hits[pos]++;
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return hits;
	}
	
	private Collection<ChipSeqLocator> collapseLocatorsByName(Collection<ChipSeqLocator> locs) { 
        LinkedHashMap<String,Map<String,Set<String>>> map = 
            new LinkedHashMap<String,Map<String,Set<String>>>();
        
        for(ChipSeqLocator loc : locs) { 
            String exptName = loc.getExptName();
            String alignName = loc.getAlignName();
            if(!map.containsKey(exptName)) { map.put(exptName, new LinkedHashMap<String,Set<String>>()); }
            if(!map.get(exptName).containsKey(alignName)) { map.get(exptName).put(alignName, new TreeSet<String>()); }
            map.get(exptName).get(alignName).addAll(loc.getReplicates());
        }
        
        LinkedList<ChipSeqLocator> collapsed = new LinkedList<ChipSeqLocator>();
        
        for(String exptName : map.keySet()) { 
            for(String alignName : map.get(exptName).keySet()) { 
                ChipSeqLocator newloc = new ChipSeqLocator(exptName, map.get(exptName).get(alignName), alignName);
                collapsed.add(newloc);
            }
        }
        
        return collapsed;
    }
	
	private RunningOverlapSum createRunningOverlapSum(Region r) {
		RunningOverlapSum totalSum = new RunningOverlapSum(r.getGenome(), r.getChrom());

		Iterator<ChipSeqHit> iter = expander.execute(r);
		while (iter.hasNext()) {
			ChipSeqHit hit = iter.next();
			Region extended = null;
			if(hit.getStrand() == '+') { 
				extended = new Region(hit.getGenome(), hit.getChrom(), hit.getStart(), hit.getEnd() + (int)readExtension);
			} else {
				extended = new Region(hit.getGenome(), hit.getChrom(), hit.getStart()-(int)readExtension, hit.getEnd());
			}
			totalSum.addRegion(extended);
		}
		return totalSum;
	}

	private StrandedRegion hit2region(int ext, ChipSeqHit hit) { 
    	if(hit.getStrand() == '+') { 
    		return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart(), hit.getEnd() + ext, hit.getStrand());
    	} else { 
    		return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart() - ext, hit.getEnd(), hit.getStrand());
    	}
    }
}
