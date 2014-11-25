package edu.psu.compbio.seqcode.deepseq.experiments;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.ReadHit;
import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.StrandedPair;
import edu.psu.compbio.seqcode.deepseq.hitloaders.*;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;

/**
 * Sample represents a single experimental sample whose hits are sourced from one or more HitLoaders. 
 * 
 * @author mahony
 *
 */
public class Sample {

	private int index;
	private Collection<HitLoader> loaders;
	private HitCache cache=null;
	private ExptConfig econfig;
	private Genome gen;
	protected String name;
	protected String sourceName=""; //String describing the source files or DBIDs 
	protected double totalHits; //totalHits is the sum of alignment weights
	protected double uniqueHits; //count of unique mapped positions (just counts the number of bases with non-zero counts - does not treat non-uniquely mapped positions differently)
	protected double totalPairs=0; //count of the total number of paired hits
	protected double uniquePairs=0; //count of the total number of unique paired hits
	protected float maxReadsPerBP=-1;
	
	
	/**
	 * Constructor
	 * @param g Genome (can be null to estimate from data)
	 * @param name String
	 */
	public Sample(int index, ExptConfig c, String name, float perBaseReadMax){
		this.index = index;
		this.name = name;
		econfig = c;
		gen=c.getGenome();
		totalHits=0;
		loaders = new ArrayList<HitLoader>();
		maxReadsPerBP= perBaseReadMax;
	}

	//Accessors
	public int getIndex(){return index;}
	public Genome getGenome(){return(gen);}
	public String getName(){return name;}
	public String getSourceName(){return sourceName;}
	public double getHitCount(){return(totalHits);}
	public double getHitPositionCount(){return(uniqueHits);}
	public double getPairCount(){return(totalPairs);}
	public double getUniquePairCount(){return(uniquePairs);}
	public void setGenome(Genome g){gen=g; cache.setGenome(g);}

	/**
	 * Add a HitLoader to the set
	 * @param h HitLoader
	 */
	public void addHitLoader(HitLoader h){
		loaders.add(h); 
		sourceName= sourceName.equals("") ? h.getSourceName() : sourceName+";"+h.getSourceName();
	}
	
	/**
	 * Initialize the cache
	 * @param cacheEntireGenome : boolean to keep the full set of hits cached
	 * @param initialCachedRegions : list of regions to keep cached at the start (can be null)
	 */
	public void initializeCache(boolean cacheEntireGenome, List<Region> initialCachedRegions){
		cache = new HitCache(econfig.getLoadPairs(), econfig, loaders, maxReadsPerBP, cacheEntireGenome, initialCachedRegions);
		totalHits = cache.getHitCount();
		uniqueHits = cache.getHitPositionCount();
		totalPairs = cache.getPairCount();
		uniquePairs = cache.getUniquePairCount();
		if(gen==null)
			gen = cache.getGenome();
	}
	
	
	/**
	 * Load all base counts in a region, regardless of strand.
	 * If caching in local files, group calls to this method by same chromosome. 
	 * @param r Region
	 * @return List of StrandedBaseCounts
	 */
	public List<StrandedBaseCount> getBases(Region r) {
		return cache.getBases(r);
	}
	/**
	 * Loads hits from a given strand in the region.
	 * If caching in local files, group calls to this method by same chromosome.
	 * @param r Region
	 * @return List of StrandedBaseCounts
	 */
	public List<StrandedBaseCount> getStrandedBases(Region r, char strand) {
		return cache.getStrandedBases(r, strand);
	}
	
	/**
	 * Load all pairs in a region, regardless of strand.
	 * If caching in local files, group calls to this method by same chromosome. 
	 * @param r Region
	 * @return List of StrandedBaseCounts
	 */
	public List<StrandedPair> getPairs(Region r) {
		return cache.getPairs(r);
	}
	
	/**
	 * Sum of all hit weights in a region.
	 * If caching in local files, group calls to this method by same chromosome.
	 * @param r Region
	 * @return float 
	 */
	public float countHits(Region r) {
		return cache.countHits(r);
	}
	/**
	 * Sum of hit weights in one strand of a region.
	 * If caching in local files, group calls to this method by same chromosome.
	 * @param r Region
	 * @return float 
	 */
    public float countStrandedBases(Region r, char strand) {
		return cache.countStrandedBases(r, strand);
    }

    /**
     * Convert all hits into ReadHits
     * @param readLen
     * @return
     */
    public List<ReadHit> exportReadHits(int readLen){
		return(cache.exportReadHits(readLen));
	}
    /**
	 * Simple count correction with a scaling factor and a floor of one. 
	 * Beware: only works if all reads are loaded.
	 * @param perBaseScaling float threshold
	 */
	public void linearCountCorrection(float perBaseScaling){
		cache.linearCountCorrection(perBaseScaling);
	}
    /**
     * Cleanup
     */
    public void close(){
    	cache.close();
    }

}
