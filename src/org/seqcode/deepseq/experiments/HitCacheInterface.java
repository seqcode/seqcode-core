package org.seqcode.deepseq.experiments;

import java.util.HashMap;
import java.util.List;

import org.seqcode.deepseq.ExtReadHit;
import org.seqcode.deepseq.ReadHit;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;

/**
 * This interface is intended to integrate the HDF5HitLoader and HDF5HitCache into the whole project
 * 
 * @author Jianyu Yang
 *
 */

public interface HitCacheInterface {
	
	//Accessors		
	public double getHitCount();
	public double getHitCountPos();
	public double getHitCountNeg();
	public double getHitPositionCount();
	public double getPairCount();
	public double getUniquePairCount();
	public Genome getGenome();
	public void setGenome(Genome g);
	
	// get the fragment size frequency distribution
	public HashMap<Integer, Integer> getFragSizeFrequency();
	
	/**
	 * Load all base counts in a region, regardless of strand.
	 * @param r Region
	 * @return List of StrandedBaseCounts
	 */
	public List<StrandedBaseCount> getBases(Region r);
	
	/**
	 * Load hits in the region on specific strand
	 * @param r
	 * @param strand
	 * @return
	 */
	public List<StrandedBaseCount> getStrandedBases(Region r, char strand);
	
	/**
	 * load pairs in the region, regardless of strand.
	 * @param r
	 * @return
	 */
	public List<StrandedPair> getPairs(Region r);
	
	/**
	 * load pairs in the region on specific strand
	 * @param r
	 * @param strand
	 * @return
	 */
	public List<StrandedPair> getPairsOnStrand(Region r, char strand);
	
	/**
	 * load pairs in the region according to the pairMid
	 * @param r
	 * @return
	 */
	public List<StrandedPair> getPairsByMid(Region r);
	
	/**
	 * count number of hits in a region
	 * @param r
	 * @return
	 */
	public float countHits(Region r);
	
	/**
	 * count number of hits in a region on specific strand
	 * @param r
	 * @param strand
	 * @return
	 */
	public float countStrandedBases(Region r, char strand);
	
	/**
	 * Convert all hits into ReadHits for a given region
	 * @param r
	 * @param readLen
	 * @return
	 */
	public List<ReadHit> exportReadHits(Region r, int readLen);
	
	/**
	 * Convert all hits into ReadHits
	 * @param readLen
	 * @return
	 */
	public List<ReadHit> exportReadHits(int readLen);
	
	/**
	 * @param r
	 * @param readLen
	 * @param startShift
	 * @param fivePrimeExt
	 * @param threePrimeExt
	 * @return
	 */
	public List<ExtReadHit> exportExtReadHits(Region r, int readLen, int startShift, int fivePrimeExt, int threePrimeExt);
	
	/**
	 * Simple count correction with a scaling factor and a floor of one.
	 * @param perBaseScaling
	 */
	public void linearCountCorrection(float perBaseScaling);
	
	/**
	 * Clean up
	 */
	public void close();
}
