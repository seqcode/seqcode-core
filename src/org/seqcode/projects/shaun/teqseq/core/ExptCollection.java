package org.seqcode.projects.shaun.teqseq.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.genome.location.Region;
import org.seqcode.gse.utils.Pair;


/**
 * ExptCollection: an object representing a collection of experimental conditions (which in turn collect sets of replicates)
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExptCollection {

	protected String name;
	protected List<ExptCondition> conditions;
	protected List<ExptReplicate> replicates;
	protected HashMap<String, Integer> replicateIndex;
	protected int rCount=0;
	protected Map<String, List<Region>> regionsOfInterest;
	
	public ExptCollection(String name, Map<String, List<Region>> regionsOfInterest){
		this.name = name;
		conditions = new ArrayList<ExptCondition>();
		replicates = new ArrayList<ExptReplicate>();
		replicateIndex = new HashMap<String, Integer>();
		this.regionsOfInterest = regionsOfInterest;
	}
	
	//Accessors
	public String getName(){return name;}
	public List<ExptCondition> getConditions(){return conditions;}
	public List<ExptReplicate> getReplicates(){return replicates;}
	public int getReplicateCount(){return rCount;}
	public int getReplicateIndex(String rName){return replicateIndex.get(rName);}
	public Map<String, List<Region>> getRegionsOfInterest(){return regionsOfInterest;}
	
	/**
	 * Add a new condition to the collection
	 * @param rep Instantiated ExptCondition
	 */
	public void addCondition(ExptCondition cond){
		conditions.add(cond);
		for(ExptReplicate r: cond.getReplicates()){
			replicateIndex.put(r.getName(), rCount);
			replicates.add(r);
			rCount++;
		}
	}
	
	/**
	 * Get AlignHit lists for the next region of interest across all experimental replicates.
	 * Assumes that all RegionHitsExtractor objects have been initialized with the same set of regions of interest.
	 * 
	 * TODO: How does this behave when a given experiment does not have any reads for one of the chromosomes?
	 * 
	 * @param seqBiasCorrected boolean: true to use the SeqBiasModel to correct for sequencing bias
	 * @return a Pair containing the current sub-chromosomal Region and a corresponding HashMap of AlignHit Lists for each experimental replicate (indexed by replicate name)
	 */
	public Pair<Region, HashMap<String, List<AlignHit>>> getAllNextRegionHit(boolean seqBiasCorrected){
		HashMap<String, List<AlignHit>> hitSet = new HashMap<String, List<AlignHit>>();
		
		Region currRegion = null;
		boolean first =true;
		for(ExptCondition c : conditions){
			for(ExptReplicate r : c.getReplicates()){
				Pair<Region, List<AlignHit>> currSet = seqBiasCorrected ? r.getNextRegionHitsWeighted() : r.getNextRegionHits();
				
				if(currSet !=null){
				    //Error check
				    if(!first){
					if(!currSet.car().equals(currRegion)){
					    System.err.println("ExptCollection: Error: All replicates did not return the same sub-chromosomal region.");
					    System.exit(1);
					}
				    }else{
					currRegion = currSet.car();
					first=false;
				    }
				    
				    //Add the current set
				    hitSet.put(r.getName(), currSet.cdr());
				}
			}
		}
		
		return(new Pair<Region, HashMap<String, List<AlignHit>>>(currRegion, hitSet));
	}
	
	/**
	 * Get AlignHit lists for the next subchromosome across all experimental replicates.
	 * Assumes that all RegionHitsExtractor objects have been initialized with the same set of regions of interest.
	 * 
	 * TODO: How does this behave when a given experiment does not have any reads for one of the chromosomes?
	 * 
	 * @param seqBiasCorrected boolean: true to use the SeqBiasModel to correct for sequencing bias
	 * @return a Pair containing the current sub-chromosomal Region and a corresponding HashMap of AlignHit Lists for each experimental replicate (indexed by replicate name)
	 */
	public Pair<Region, HashMap<String, List<AlignHit>>> getAllNextSubChrHit(boolean seqBiasCorrected){
		HashMap<String, List<AlignHit>> hitSet = new HashMap<String, List<AlignHit>>();
		
		Region currRegion = null;
		boolean first =true;
		for(ExptCondition c : conditions){
			for(ExptReplicate r : c.getReplicates()){
				Pair<Region, List<AlignHit>> currSet = seqBiasCorrected ? r.getNextSubChrHitsWeighted() : r.getNextSubChrHits();
				
				if(currSet !=null){
				    //Error check
				    if(!first){
						if(!currSet.car().equals(currRegion)){
						    System.err.println("ExptCollection: Error: All replicates did not return the same sub-chromosomal region.");
						    System.exit(1);
						}
				    }else{
						currRegion = currSet.car();
						first=false;
				    }
				    
				    //Add the current set
				    hitSet.put(r.getName(), currSet.cdr());
				}
			}
		}
		
		return(new Pair<Region, HashMap<String, List<AlignHit>>>(currRegion, hitSet));
	}
	
	/**
	 * Checks that all experimental replicates have regions remaining in the RegionHitsExtractor
	 * TODO: How does this behave when a given experiment does not have any reads for one of the chromosomes?
	 * @return boolean
	 */
	public boolean hasNextRegions(){
		boolean hasNext = true;
		for(ExptCondition c : conditions){
			for(ExptReplicate r : c.getReplicates()){
				hasNext = hasNext & r.hasNextRegion();
			}
		}
		return(hasNext);
	}
	
	/**
	 * Reset the read loaders
	 */
	public void resetHitExtractors(){
		for(ExptCondition c : conditions){
			c.resetHitExtractors();
		}
	}
}
