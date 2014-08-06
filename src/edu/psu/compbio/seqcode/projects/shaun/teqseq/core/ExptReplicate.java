package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.ObsExpSeqBiasModel;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.SeqBiasModel;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AGene;

/**
 * ExptReplicate: an object representing a single replicate of a given experimental condition. 
 * Contains a read loader and various objects required to normalize a count or correct biases. 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExptReplicate {

	protected String name;
	protected int ID;
	protected GenomeLoader gLoad;
	protected ReadLoader reads;
	protected RegionHitsExtractor hitsExtractor=null;
	protected Map<String, List<Region>> regionsOfInterest;
	protected SeqBiasModel seqBias;
	protected double scalingFactor=1.0;
	
	/**
	 * Constructor
	 * @param gl GenomeLoader
	 * @param rl ReadLoader
	 * @param initialRegions Collections of Regions of interest organized by chromosome name
	 * @param genes Set of gene annotations only used by the SeqBiasModel
	 */
	public ExptReplicate(String name, int id, GenomeLoader gl, ReadLoader rl, Map<String, List<Region>> regionsOfInterest, Collection<AGene> genes, boolean correctSeqBias){
		this.name = name;
		this.ID = id;
		gLoad = gl;
		reads = rl;
		
		System.err.println("Adding experiment replicate "+name);
		if(correctSeqBias){
			seqBias = new ObsExpSeqBiasModel(gLoad, 6, regionsOfInterest);
			if(!seqBias.attemptLoad(reads)){
				if(gLoad.seqIsAvailable()){
					System.err.println("ExptReplicate: Cannot find sequencing bias model for "+name+", training a new model.");
					seqBias.execute(reads, genes);
				}else{
					System.err.println("ExptReplicate: Cannot find sequencing bias model for "+name+" and genome sequence was not provided for training the model.");
					seqBias=null;
				}
			}	
		}else{
			seqBias=null;
		}
		initializeHitExtractor(regionsOfInterest);
	}
	
	//Accessors
	public String getName(){return name;}
	public int getID(){return ID;}
	
	/**
	 * Initialize or update the RegionHitsExtractor object with a set of regions of interest
	 * @param regions Collections of Regions of interest organized by chromosome name
	 */
	public void initializeHitExtractor(Map<String, List<Region>> regions){
		regionsOfInterest = regions;
		if(hitsExtractor!=null)
			hitsExtractor.close();
		hitsExtractor = new RegionHitsExtractor(gLoad, reads, regionsOfInterest,regionsOfInterest);
	}
	
	/**
	 * Get a set of weighted hits for a region of interest from the hits extractor
	 * @return a Pair containing a Region and associated hits
	 */
	public Pair<Region, List<AlignHit>> getNextRegionHitsWeighted(){
		return hitsExtractor.getNext(seqBias);
	}
	
	/**
	 * Get a set of unweighted hits for a region of interest from the hits extractor
	 * @return a Pair containing a Region and associated hits
	 */
	public Pair<Region, List<AlignHit>> getNextRegionHits(){
		return hitsExtractor.getNext();
	}
	
	/**
	 * Get a set of weighted hits for a sub-chromosome from the hits extractor
	 * @return a Pair containing a Region and associated hits
	 */
	public Pair<Region, List<AlignHit>> getNextSubChrHitsWeighted(){
		return hitsExtractor.getNextSubChr(seqBias);
	}
	
	/**
	 * Get a set of unweighted hits for a sub-chromosome from the hits extractor
	 * @return a Pair containing a Region and associated hits
	 */
	public Pair<Region, List<AlignHit>> getNextSubChrHits(){
		return hitsExtractor.getNextSubChr();
	}
	
	/**
	 * Query the hits extractor to see if there are any reads remaining in the read iterator. 
	 * @return boolean
	 */
	public boolean hasNextRegion(){
		return hitsExtractor.hasNextRegion();
	}
	
	/**
	 * Reset the hit extractor
	 */
	public void resetHitExtractor(){
		hitsExtractor.resetReadLoader();
	}
}
