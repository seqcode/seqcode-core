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

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.RunningOverlapSum;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqData {
    public static final int defaultReadLength = 26, defaultReadExtension = 174;
	private double readLength = defaultReadLength;
	private double readExtension = defaultReadExtension;
	private String exptName;
	private Genome currentGen=null;
	private double hitCount = 0;
	private double totalSeq = 0;
	private SeqExpander expander=null;
	private SeqLocator loc = null;
	private SeqDataLoader loader;
	private LinkedList<SeqAlignment> alignments;

	public SeqData(Genome g, String exptName){
		this(g, exptName, "");
	}
	 
	public SeqData(Genome g,  String exptName, String replicate){
		currentGen=g;
		this.exptName=exptName;
		
		SeqDataLoader chipSeqLoader = null;
    	try {
    		chipSeqLoader = new SeqDataLoader();
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
        			List<SeqLocator> locs = new Vector<SeqLocator>();
        			
        			List<SeqExpt> expts = new Vector<SeqExpt>();
        			expts.addAll(chipSeqLoader.loadExperiments(exptName));
					
        			for (SeqExpt expt : expts) {
        				Collection<SeqAlignment> aligns;
						aligns = chipSeqLoader.loadAllAlignments(expt);
						
                		for (SeqAlignment currentAlign : aligns) {
                			if (currentAlign.getGenome().equals(g)) { 
                				SeqLocator currentLoc = new SeqLocator(expt.getName(), 
                                        expt.getReplicate(), currentAlign.getName());
                				locs.add(currentLoc);
        						break;
        					}
                		}
        			}

        			List<SeqLocator> collapsedLocs = new Vector<SeqLocator>(this.collapseLocatorsByName(locs));
        			if (collapsedLocs.size() != 1) {
        				System.err.println(collapsedLocs.size() + " collapsed locators");
        				System.exit(0);
        			}
        			loc = collapsedLocs.get(0);
        		}
        		else {
        			SeqExpt expt;
					expt = chipSeqLoader.loadExperiment(exptName, replicate);
					
        			SeqAlignment align = null;
            		Collection<SeqAlignment> aligns;
					aligns = chipSeqLoader.loadAllAlignments(expt);
					
            		for (SeqAlignment currentAlign : aligns) {
            			if (currentAlign.getGenome().equals(g)) { 
    						align = currentAlign;
    						break;
    					}
            		}
            		
            		loc = new SeqLocator(expt.getName(), expt.getReplicate(), align.getName());
        		}
        		
        		expander = new SeqExpander(loc);
        	    loader = new SeqDataLoader();
    	        alignments = new LinkedList<SeqAlignment>();
    	        
    	        if(loc.getReplicates().isEmpty()) { 
    	        	Collection<SeqExpt> expts = loader.loadExperiments(loc.getExptName());
    	        	for(SeqExpt expt : expts) { 
    	        		SeqAlignment alignment = 
    	        			loader.loadAlignment(expt, loc.getAlignName(), g);
    	        		if(alignment != null) { 
    	        			alignments.add(alignment);
    	        		}
    	        	}
    	        } else {
    	            for(String repName : loc.getReplicates()) { 
    	                SeqExpt expt = loader.loadExperiment(loc.getExptName(), repName);
    	                SeqAlignment alignment = 
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
    public SeqData(Genome g, SeqLocator locator) throws NotFoundException, SQLException {
		currentGen=g;
        exptName = locator.getExptName();
        try {
			loader = new SeqDataLoader();
			
	        alignments = new LinkedList<SeqAlignment>();
	        if (locator.getAlignName() == null) {
	            if(locator.getReplicates().isEmpty()) { //No alignment name, no replicate names
	            	Collection<SeqExpt> expts = loader.loadExperiments(locator.getExptName());
	        		List<SeqLocator> locs = new Vector<SeqLocator>();
	        		for(SeqExpt expt : expts) { 
	                	Collection<SeqAlignment> aligns;
						aligns = loader.loadAllAlignments(expt);
						for (SeqAlignment currentAlign : aligns) {
	            			if (currentAlign.getGenome().equals(g)) { 
	            				SeqLocator currentLoc = new SeqLocator(expt.getName(), 
	                                    expt.getReplicate(), currentAlign.getName());
	            				locs.add(currentLoc);
	            				alignments.add(currentAlign);
	    						break;
	    					}
	            		}
	    			}
	    			List<SeqLocator> collapsedLocs = new Vector<SeqLocator>(this.collapseLocatorsByName(locs));
	    			if (collapsedLocs.size() != 1) {
	    				System.err.println(collapsedLocs.size() + " collapsed locators");
	    				System.exit(0);
	    			}
	    			locator = collapsedLocs.get(0);
	            } else { //No alignment name, given replicate names
	                for(String repName : locator.getReplicates()) { 
	                    SeqExpt expt = loader.loadExperiment(locator.getExptName(), repName);
	                    SeqAlignment alignment = 
	                        loader.loadAlignment(expt, locator.getAlignName(), g);
	                    if(alignment != null) { 
	                        locator = new SeqLocator(locator.getExptName(),
	                                                     locator.getReplicates(),
	                                                     alignment.getName());
	                        alignments.add(alignment);
	                        break;
	                    }
	                }
	            }
	        } else {
	        	if(locator.getReplicates().isEmpty()) {//Given alignment name, no replicate names
	        		Collection<SeqExpt> expts = loader.loadExperiments(locator.getExptName());
	        		List<SeqLocator> locs = new Vector<SeqLocator>();
	        		for(SeqExpt expt : expts) { 
	                	Collection<SeqAlignment> aligns;
						aligns = loader.loadAllAlignments(expt);
						for (SeqAlignment currentAlign : aligns) {
	            			if (currentAlign.getGenome().equals(g) && currentAlign.getName().equals(locator.getAlignName())) { 
	            				SeqLocator currentLoc = new SeqLocator(expt.getName(), 
	                                    expt.getReplicate(), currentAlign.getName());
	            				locs.add(currentLoc);
	            				alignments.add(currentAlign);
	    						break;
	    					}
	            		}
	    			}
	    			List<SeqLocator> collapsedLocs = new Vector<SeqLocator>(this.collapseLocatorsByName(locs));
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


	        expander = new SeqExpander(locator);
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
	public SeqExpander getExpander(){return expander;}
	public SeqLocator getLocator(){return loc;}
	
	
	public double countHits(){
		hitCount=0;
		totalSeq=0;
		Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			
			hitCount += expander.getHitCount(currentChrom);
			totalSeq += (double)currentChrom.getWidth();
		}		
		
		return hitCount;
	}
	public double countHits(Region a) {
		double hits=0;
		try {
			for(SeqAlignment alignment : alignments) { 
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
			for(SeqAlignment alignment : alignments) { 
				Collection<SeqHit> hits = loader.loadByRegion(alignment, a);
				for(SeqHit curr : hits){
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
			for(SeqAlignment alignment : alignments) { 
				Collection<SeqHit> hits = loader.loadByRegion(alignment, a);
				for(SeqHit curr : hits){
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
			for(SeqAlignment alignment : alignments) { 
				Iterator<SeqHit> h = loader.loadByRegion(alignment, a).iterator();
				while(h.hasNext()){
					SeqHit csh = h.next();
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
	
	private Collection<SeqLocator> collapseLocatorsByName(Collection<SeqLocator> locs) { 
        LinkedHashMap<String,Map<String,Set<String>>> map = 
            new LinkedHashMap<String,Map<String,Set<String>>>();
        
        for(SeqLocator loc : locs) { 
            String exptName = loc.getExptName();
            String alignName = loc.getAlignName();
            if(!map.containsKey(exptName)) { map.put(exptName, new LinkedHashMap<String,Set<String>>()); }
            if(!map.get(exptName).containsKey(alignName)) { map.get(exptName).put(alignName, new TreeSet<String>()); }
            map.get(exptName).get(alignName).addAll(loc.getReplicates());
        }
        
        LinkedList<SeqLocator> collapsed = new LinkedList<SeqLocator>();
        
        for(String exptName : map.keySet()) { 
            for(String alignName : map.get(exptName).keySet()) { 
                SeqLocator newloc = new SeqLocator(exptName, map.get(exptName).get(alignName), alignName);
                collapsed.add(newloc);
            }
        }
        
        return collapsed;
    }
	
	private RunningOverlapSum createRunningOverlapSum(Region r) {
		RunningOverlapSum totalSum = new RunningOverlapSum(r.getGenome(), r.getChrom());

		Iterator<SeqHit> iter = expander.execute(r);
		while (iter.hasNext()) {
			SeqHit hit = iter.next();
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

	private StrandedRegion hit2region(int ext, SeqHit hit) { 
    	if(hit.getStrand() == '+') { 
    		return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart(), hit.getEnd() + ext, hit.getStrand());
    	} else { 
    		return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart() - ext, hit.getEnd(), hit.getStrand());
    	}
    }
}
