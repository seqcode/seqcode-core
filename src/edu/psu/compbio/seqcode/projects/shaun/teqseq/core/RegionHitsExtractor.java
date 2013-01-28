package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.SeqBiasModel;

/**
 * Extract read hits for a set of regions of interest. We assume that the read stream is pre-sorted. 
 * We also assume that the regions of interest are also sorted and organized by chromosome. 
 * This class splits the chromosomes up into large non-overlapping sub-regions. 
 * For each sub-region, all corresponding AlignHits are loaded to memory from the ReadLoader. 
 * Then, for each region of interest, corresponding sub-lists of AlignHits are loaded. 
 * This approach allows us to load reads from the stream while querying regions that can overlap each other.
 * The sub-chromosomal regions are concatenations of the regions of interest - 
 * they are at least 5Mbp wide, and at least 500Kbp of space has to exist between neighboring regions of interest. 
 *   (Known Ensembl genes are up to 4.5Mbp wide with up to 600Kbp between exons)
 * 
 * This class can be used to iterate over the regions of interest, loading each set of corresponding align hits in turn.
 * Alternatively, you can iterate directly over the sub-chromosomal regions, returning all align hits in that large region.
 * I can't guarantee what will happen if you try both types of queries in the same iteration. Nothing good will come of it anyway. 
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class RegionHitsExtractor {

	protected GenomeLoader gLoad;
	protected ReadLoader reads;
	protected Map<String, List<Region>> regionsOfInterest;
	protected Map<String, List<Region>> subChromosomes;
	protected int currSubChrNum=0;
	protected List<AlignHit> subHitList=null;
	protected String currChr="";
	protected int currRegNum;
	protected boolean subListInitialized=false;
	
	/**
	 * Constructor
	 * @param g GenomeLoader
	 * @param rl ReadLoader
	 * @param cRegs Collections of Regions of interest organized by chromosome name
	 */
	public RegionHitsExtractor(GenomeLoader g, ReadLoader rl, Map<String, List<Region>> interestRegs, Map<String, List<Region>> subChrRegs){
		gLoad = g;
		reads = rl;
		regionsOfInterest = interestRegs;
		subChromosomes = subChrRegs;
		currRegNum=0; currSubChrNum=0;
		reads.initializeReadIterator();
		if(reads.iteratorHasNextHit()){
			currChr=reads.getCurrHit().getChrom();
			currRegNum=0;
		}
		
		/*
		int numRegInterest=0, numSubChr=0;
		for(String s : regionsOfInterest.keySet())
			numRegInterest+=regionsOfInterest.get(s).size();
		for(String s : subChromosomes.keySet())
			numSubChr+=subChromosomes.get(s).size();
		System.err.println("RegionHitsExtractor: Initialized with "+numSubChr+" sub-chromosomal regions.");
		*/
	}

	/**
	 * Get the next pair of Regions and corresponding List of AlignHits
	 * @param bias SeqBiasModel (trained) that gives the relative weights of each sequence. This can be null.
	 * @return a Pair containing a Region and a corresponding List of AlignHits	 * 
	 */
	public Pair<Region, List<AlignHit>> getNext(SeqBiasModel bias){
		if(!subListInitialized){
			//Initialize the hit list for this chromosome if there are regions of interest therein
			if(subChromosomes.containsKey(currChr)){
				subHitList = reads.getOverlappingHitsFromStream(subChromosomes.get(currChr).get(currSubChrNum), bias);
				//System.err.println("Initialized with sub-chromosome "+subChromosomes.get(currChr).get(currSubChrNum).getLocationString()+", which has "+subHitList.size()+" hits.");
			}
			subListInitialized=true;			
		}
		
		//currChr is determined based on the position of the read iterator
		//If there is no region of interest in the current chromosome, 
		// or if we have gone past the end of a chromosome, or processed all regions of interest
		//  fast-forward the read iterator
		while((!subChromosomes.containsKey(currChr) || currSubChrNum>=subChromosomes.get(currChr).size() || currRegNum>=regionsOfInterest.get(currChr).size()) && reads.iteratorHasNextHit()){
			while(currChr.equals(reads.getCurrHit().getChrom()) && reads.nextRead()){}
			if(reads.iteratorHasNextHit()){
				currChr=reads.getCurrHit().getChrom();
				currRegNum=0;
				currSubChrNum=0;
				if(subChromosomes.containsKey(currChr)){
					//load sublist
					subHitList = reads.getOverlappingHitsFromStream(subChromosomes.get(currChr).get(currSubChrNum), bias);
					//System.err.println("Loaded sub-chromosome "+subChromosomes.get(currChr).get(currSubChrNum).getLocationString()+", which has "+subHitList.size()+" hits.");
				}
			}
		}
		if(subChromosomes.containsKey(currChr)){
			
			//If there is another region of interest 
			if(currSubChrNum<subChromosomes.get(currChr).size() && currRegNum<regionsOfInterest.get(currChr).size()){
				//currRegion is the next region for which we need reads
				Region currRegion = regionsOfInterest.get(currChr).get(currRegNum);
			
				//If the current region of interest is not in the current sub-chromosome 
				if(!subChromosomes.get(currChr).get(currSubChrNum).contains(currRegion)){
					//Load the next sub-chromosome
					currSubChrNum++;
					if(currSubChrNum<subChromosomes.get(currChr).size()){
						//Load the current sub-chromosome's reads
						subHitList = reads.getOverlappingHitsFromStream(subChromosomes.get(currChr).get(currSubChrNum), bias);
						//System.err.println("Loaded sub-chromosome "+subChromosomes.get(currChr).get(currSubChrNum).getLocationString()+", which has "+subHitList.size()+" hits.");
					}else{
						return null;
					}
				}
				//System.err.println("RegionHitsExtractor: Loading region "+currRegion.getLocationString()+" from sub-chromosome "+subChromosomes.get(currChr).get(currSubChrNum).getLocationString()+", which has "+subHitList.size()+" hits.");
				//Finally, get the reads
				List<AlignHit> hits = getHitsFromSubList(currRegion, subHitList);
				currRegNum++;
				System.err.println("RegionHitsExtractor: Loaded region "+currRegion.getLocationString()+" with "+hits.size()+" alignment hits.");
				return(new Pair<Region, List<AlignHit>>(currRegion, hits));
			}else{
				return null;
			}
		}else{
			return null;
		}
	}
	public Pair<Region, List<AlignHit>> getNext(){return this.getNext(null);}
	
	
	/**
	 * Get the next pair of Sub-Chromosomal Region and corresponding List of AlignHits
	 * Note that in this method, we ignore the regions of interest and just return all AlignHits in a sub-chromosome
	 * @param bias SeqBiasModel (trained) that gives the relative weights of each sequence. This can be null.
	 * @return a Pair containing a Region and a corresponding List of AlignHits	 * 
	 */
	public Pair<Region, List<AlignHit>> getNextSubChr(SeqBiasModel bias){
		//currChr is determined based on the position of the read iterator
		//Ignore regions of interest and focus on the sub-chromosomes
		while((!subChromosomes.containsKey(currChr) || currSubChrNum>=subChromosomes.get(currChr).size()) && reads.iteratorHasNextHit()){
			while(currChr.equals(reads.getCurrHit().getChrom()) && reads.nextRead()){}
			if(reads.iteratorHasNextHit()){
				currChr=reads.getCurrHit().getChrom();
				currRegNum=0;
				currSubChrNum=0;
			}
		}
		if(subChromosomes.containsKey(currChr) && reads.iteratorHasNextHit()){
			if(currSubChrNum<subChromosomes.get(currChr).size()){
				subHitList = reads.getOverlappingHitsFromStream(subChromosomes.get(currChr).get(currSubChrNum), bias);
				Pair<Region, List<AlignHit>> p = new Pair<Region, List<AlignHit>>(subChromosomes.get(currChr).get(currSubChrNum), subHitList);
				currSubChrNum++;
				return(p);
			}else{
			    return null;
			}
		}else{			    
		    return null;
		}
	}
	public Pair<Region, List<AlignHit>> getNextSubChr(){return this.getNextSubChr(null);}
	
	/**
	 * 
	 * @return boolean true if the read iterator has hits left
	 */
	public boolean hasNextRegion(){
	    if(reads.iteratorHasNextHit()){// && ((subChromosomes.containsKey(currChr) && currSubChrNum<subChromosomes.get(currChr).size()) && (regionsOfInterest.containsKey(currChr) && currRegNum<regionsOfInterest.get(currChr).size())))
		return true;
	    }
	    return false;
	}
	
	/**
	 * Get hits from the input list that overlap a region
	 * @param currReg the Region of interest
	 * @param hits AlignHits to filter
	 * @return List of AlignHits
	 */
	private List<AlignHit> getHitsFromSubList(Region currReg, List<AlignHit> hits){
		List<AlignHit> currHits = new ArrayList<AlignHit>();
		int i = Math.abs(Collections.binarySearch(hits, currReg))-1; //Binarysearch should return (-(insertion point)-1) when the item is not in the list
		while(i<hits.size() && currReg.contains(hits.get(i))){
			currHits.add(hits.get(i));
			i++;
		}
		return(currHits);
	}
	
	/**
	 * Reset the read iterator for another cycle through the data
	 */
	public void resetReadLoader(){
	    currRegNum=0; currSubChrNum=0;
	    reads.initializeReadIterator();
	    if(reads.iteratorHasNextHit()){
		currChr=reads.getCurrHit().getChrom();
		currRegNum=0;
	    }
	    subListInitialized=false;
	}
	
	/**
	 * Close the region hit extractor
	 */
	public void close(){
		reads.closeReader();
	}
}
