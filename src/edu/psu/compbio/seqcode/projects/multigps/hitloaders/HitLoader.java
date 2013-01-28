package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import edu.psu.compbio.seqcode.projects.multigps.framework.Read;
import edu.psu.compbio.seqcode.projects.multigps.framework.ReadHit;

/**
 * HitLoaders load alignment hits from various sources, including ReadDB and various files.
 * Five-prime positions and associated weight sums are loaded into ArrayLists. 
 * Where/how those hits are sourced is implementation-specific. 
 * This class combines functionality from ReadLoaders, AlignmentFileReaders, and ReadCache in the old setup.  
 * 
 * Five prime positions and weights are loaded into two collections of ArrayLists, where the collections are indexed by chromosome name. 
 * Within each chromosome's set, a 2D array of ArrayLists collects data for each strand.
 * However, the ArrayLists are temporary -- once a Sample loads the hits into primitive arrays, the ArrayLists are reset and the 
 * garbage collector is called. 
 * 
 * @author shaun
 *
 */
public abstract class HitLoader {

	protected double totalHits; //totalHits is the sum of alignment weights
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
	 */
	private HashMap<String, ArrayList<Float>[]> fivePrimeCountsList = null;
		
	/**
	 * Constructor
	 * @param g Genome
	 */
	public HitLoader(){
		totalHits=0;		
	}

//Accessors
	public double getHitCount(){return(totalHits);}
	public HashMap<String, ArrayList<Integer>[]> getFivePrimePositions(){return fivePrimePosList;}
	public HashMap<String, ArrayList<Float>[]> getFivePrimeCounts(){return fivePrimeCountsList;}
	
//Abstract methods
	/**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 */
	public abstract void sourceReads();

	
//Shared methods
	/**
	 * Initialize the genome and data structures. Source hits for the lists
	 */
	public void initialize(){
		resetLoader();
		
		fivePrimePosList = new HashMap<String, ArrayList<Integer>[]>();
		fivePrimeCountsList = new HashMap<String, ArrayList<Float>[]>();
	}
	
	/**
	 * Reset the loaders -- destroy the lists and call the garbage collector
	 */
	public void resetLoader(){
		fivePrimePosList=null;
		fivePrimeCountsList=null; 
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
	}
	
	/**
	 * Perform any necessary cleanup. For ReadDB, this means close the clients.
	 */
	public abstract void cleanup();
}
