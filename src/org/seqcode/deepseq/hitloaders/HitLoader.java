package org.seqcode.deepseq.hitloaders;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import org.seqcode.deepseq.HitPair;
import org.seqcode.deepseq.Read;
import org.seqcode.deepseq.ReadHit;


/**
 * HitLoaders load alignment hits & pairs from various sources, including ReadDB and various files.
 * Five-prime positions and associated weight sums are loaded into ArrayLists. 
 * Pairing information is loaded if requested and if it exits. 
 * Where/how those hits & pairs are sourced is implementation-specific. 
 * 
 * Five prime positions and weights are loaded into two collections of ArrayLists, where the collections are indexed by chromosome name. 
 * Within each chromosome's set, a 2D array of ArrayLists collects data for each strand.
 * However, the ArrayLists are temporary -- once a Sample loads the hits into primitive arrays, the ArrayLists are reset and the 
 * garbage collector is called. 
 * 
 * @author mahony
 * This class combines functionality from ReadLoaders, AlignmentFileReaders, and ReadCache in the old setup.
 */
public abstract class HitLoader {

	protected boolean loadType1=true; //Load type1 reads
	protected boolean loadType2=false; //Load type2 reads (if exists)
	protected boolean loadRead2=true; //Load read 2 in paired-end
	protected boolean loadPairs=false; //Load pair information (if exists)
	protected boolean hasPairs = false; //Flag to say there are pairs in the sample 
	protected double totalHits; //totalHits is the sum of alignment weights
	protected String sourceName=""; //String describing the source
	
	/**
	 * Five prime ends of the read hits. <br>
	 * HashMap is indexed by chromosome name. <br>
	 * Dimension in the array of ArrayLists represents the strand. 0 for '+', 1 for '-' 
	 */
	private HashMap<String, ArrayList<Integer>[]> fivePrimePosList = null;
	/**
	 * Sum of read hit weights that corresponds to the 5' position
	 * HashMap is indexed by chromosome name. <br>
	 * Dimension in the array of ArrayLists represents the strand. 0 for '+', 1 for '-'
	 * Ordering of each ArrayList is the same as fivePrimePosList 
	 */
	private HashMap<String, ArrayList<Float>[]> fivePrimeCountsList = null;
	/**
	 * R2 read hit pairing information for each R1 read hit (if pairs exist)
	 * HashMap is indexed by R1 read chromosome name. <br>
	 * Dimension in the array of ArrayLists represents the R1 read strand. 0 for '+', 1 for '-'
	 * Ordering of each ArrayList is the same as fivePrimePosList.   
	 * 
	 */
	private HashMap<String, ArrayList<HitPair>[]> hitPairsList = null;
	
		
	/**
	 * Constructor
	 * @param g Genome
	 */
	public HitLoader(boolean loadT1, boolean loadT2, boolean loadRead2, boolean loadPairs){
		this.loadType1=loadT1;
		this.loadType2=loadT2;
		this.loadRead2 = loadRead2;
		this.loadPairs=loadPairs;
		totalHits=0;		
	}

//Accessors
	public String getClassName() {return this.getClass().getSimpleName();}
	public boolean hasPairedReads(){return hasPairs;}
	public double getHitCount(){return(totalHits);}
	public String getSourceName(){return sourceName;}
	public HashMap<String, ArrayList<Integer>[]> getFivePrimePositions(){return fivePrimePosList;}
	public HashMap<String, ArrayList<Float>[]> getFivePrimeCounts(){return fivePrimeCountsList;}
	public HashMap<String, ArrayList<HitPair>[]> getPairs(){return hitPairsList;}
	
//Abstract methods
	/**
	 * Get all hits from the appropriate source (implementation-specific).
	 * Loads single end data to the fivePrimePosList and fivePrimeCountsList.
	 * Enforcing which reads to load (Type1 and/or Type2) is also implementation-specific. 
	 * Loads pairs to hitPairsList (if requested & if they exist).
	 * 
	 */
	public abstract void sourceAllHits();

	
//Shared methods
	/**
	 * Initialize the genome and data structures. Source hits for the lists
	 */
	public void initialize(){
		resetLoader();
		
		fivePrimePosList = new HashMap<String, ArrayList<Integer>[]>();
		fivePrimeCountsList = new HashMap<String, ArrayList<Float>[]>();
		if(loadPairs)
			hitPairsList = new HashMap<String, ArrayList<HitPair>[]>();
	}
	
	/**
	 * Reset the loaders -- destroy the lists and call the garbage collector
	 */
	public void resetLoader(){
		//Free memory
		if(fivePrimePosList!=null){
			for(String chr: fivePrimePosList.keySet()){
				fivePrimePosList.get(chr)[0].clear();
				fivePrimePosList.get(chr)[1].clear();
			}
			fivePrimePosList.clear();
		}
		if(fivePrimeCountsList!=null){
			for(String chr: fivePrimeCountsList.keySet()){
				fivePrimeCountsList.get(chr)[0].clear();
				fivePrimeCountsList.get(chr)[1].clear();
			}
			fivePrimeCountsList.clear();
		}
		if(loadPairs && hitPairsList!=null){
			for(String chr: hitPairsList.keySet()){
				hitPairsList.get(chr)[0].clear();
				hitPairsList.get(chr)[1].clear();
			}
			hitPairsList.clear();
		}
		System.gc();
	}
	
	/**
	 * 	Add hits to the list data structures.
	 * 	It may be called multiple times to retrieve all the data, then populateArrays() is called 
	 */
	protected void addHits(String chrom, char strand, Collection<Integer> coords, Collection<Float> counts){
		int strandInd = strand == '+' ? 0 : 1;
		if(!fivePrimePosList.containsKey(chrom))
			addChr(chrom);
		fivePrimePosList.get(chrom)[strandInd].addAll(coords);
		fivePrimeCountsList.get(chrom)[strandInd].addAll(counts);
		for (float c: counts)
			totalHits += c;
	}//end of addHits method
	
	/**
	 * Add hits to the list data structures from a Read
	 * @param r Read
	 */
	protected void addHits(Read r){
		for(ReadHit h : r.getHits()){
			char strand = h.getStrand();
			int strandInd = strand == '+' ? 0 : 1;
			if(!fivePrimePosList.containsKey(h.getChrom()))
				addChr(h.getChrom());
			fivePrimePosList.get(h.getChrom())[strandInd].add(strand == '+' ?h.getStart():h.getEnd());
			fivePrimeCountsList.get(h.getChrom())[strandInd].add(h.getWeight());
			totalHits++;
		}
	}//end of addHits method
	
	/**
	 * Add paired hit information to the list data structure
	 * @param HitPair collection
	 */
	protected void addPairs(String chrom, char strand, Collection<HitPair> pairs){
		if(!hasPairs){
			//This is the first pair being added.
			hasPairs=true;
		}
		int strandInd = strand == '+' ? 0 : 1;
		if(!hitPairsList.containsKey(chrom))
			addChr(chrom);
		hitPairsList.get(chrom)[strandInd].addAll(pairs);
	}
	/**
	 * Add paired hit information to the list data structure
	 * @param HitPair
	 */
	protected void addPair(String chrom, char strand, HitPair pair){
		if(!hasPairs){
			//This is the first pair being added.
			hasPairs=true;
		}
		int strandInd = strand == '+' ? 0 : 1;
		if(!hitPairsList.containsKey(chrom))
			addChr(chrom);
		hitPairsList.get(chrom)[strandInd].add(pair);
	}
	
	/**
	 * Add a chromosome to the hit lists
	 * @param chr String
	 */
	protected void addChr(String chr){
		ArrayList<Integer>[] currIArrayList = new ArrayList[2];
		currIArrayList[0]=new ArrayList<Integer>();
		currIArrayList[1]=new ArrayList<Integer>();
		fivePrimePosList.put(chr, currIArrayList);
		ArrayList<Float>[] currFArrayList = new ArrayList[2];
		currFArrayList[0]=new ArrayList<Float>();
		currFArrayList[1]=new ArrayList<Float>();
		fivePrimeCountsList.put(chr, currFArrayList);
		if(loadPairs){
			ArrayList<HitPair>[] currPArrayList = new ArrayList[2];
			currPArrayList[0]=new ArrayList<HitPair>();
			currPArrayList[1]=new ArrayList<HitPair>();
			hitPairsList.put(chr, currPArrayList);
		}
	}
	
	/**
	 * Perform any necessary cleanup. For ReadDB, this means close the clients.
	 */
	public abstract void cleanup();
	
}
