package org.seqcode.deepseq.experiments;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.optim.nonlinear.vector.Weight;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.hitloaders.HDF5HitLoader;
import org.seqcode.deepseq.utils.HierarchicalHitInfo;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;

public class HDF5HitCache implements HitCacheInterface{
	
	private HierarchicalHitInfo readHHI;
	private HierarchicalHitInfo pairHHI;
	private Collection<HDF5HitLoader> loaders;	// source of hits
	private File localCacheDir;
	private Genome gen;
	private ExptConfig econfig;
	protected double totalHits = 0;		// totalHits is the sum of alignment weights
	protected double totalHitsPos = 0;	// totalHitsPos is the sum of alignment weights on the plus strand
	protected double totalHitsNeg = 0;	// totalHitsNeg is the sum of alignment weights on the minus strand
	protected double uniqueHits = 0;	// count of unique mapped positions (just counts the number of bases with non-zero counts - does not treat non-uniquely mapped positions differently)
	protected double totalPairs = 0;	// count of the total number of paired hits
	protected double uniquePairs = 0;	// count of the total number of unique paired hits
	protected boolean loadReads;		// load single read
	protected boolean loadPairs;		// load them if they exist
	protected boolean hasPairs;			// some pairs have been loaded
	
	protected String readReference = "5'end";	
	protected String pairReference = "pairMid";
	
	public HDF5HitCache(boolean loadReads, boolean loadPairs, ExptConfig ec, Collection<HDF5HitLoader> hloaders) {
		this.econfig = ec;
		this.gen = econfig.getGenome();
		this.loaders = hloaders;
		
		this.loadReads = loadReads;
		this.loadPairs = loadPairs;
		
		this.localCacheDir = new File("./");
		
		initialize();
	}
	
	/**
	 * load hit info of all hitloaders to a single dataset and sort them by reference element
	 */
	private void initialize() {
		
		// Create the hitInfo dataset for the HDF5HitCache itself
		readHHI = new HierarchicalHitInfo(gen, localCacheDir.getAbsolutePath()+"HDF5HitCache.read.h5", false);
		pairHHI = new HierarchicalHitInfo(gen, localCacheDir.getAbsolutePath()+"HDF5HitCache.pair.h5", true);
		
		// Initialize
		readHHI.initializeHDF5();
		pairHHI.initializeHDF5();
		
		// Append all hits into the dataset and extract the info into the HDF5HitCache
		for(HDF5HitLoader currLoader: loaders) {
			currLoader.sourceAllHits();
			
			// readHHI
			HierarchicalHitInfo currloaderReadHHI = currLoader.getReadInfo();
			for (String chrom: gen.getChromList())
				for (int strand=0; strand<2; strand++){
					int startOrder = 0;
					int length = currloaderReadHHI.getLength(chrom, strand);
					for(int endOrder=10000; endOrder<=length; endOrder+=10000) {
						try {
							readHHI.appendArray(chrom, strand, currloaderReadHHI.getElement(chrom, strand, startOrder, endOrder));
						} catch (Exception e) {
							e.printStackTrace();
						}
						startOrder = endOrder;
					}
					// append the rest part of the hits into the hitcache
					if (currloaderReadHHI.getLength(chrom, strand)>startOrder)
						try {
							readHHI.appendArray(chrom, strand, currloaderReadHHI.getElement(chrom, strand, startOrder, length));
						} catch (Exception e) {
							e.printStackTrace();
						}
				}
			
			// pairHHI
			HierarchicalHitInfo currloaderPairHHI = currLoader.getHitPairInfo();
			for (String chrom: gen.getChromList())
				for (int strand=0; strand<2; strand++){
					int startOrder = 0;
					int length = currloaderPairHHI.getLength(chrom, strand);
					if(length > 0) {hasPairs = true;}
					for(int endOrder=10000; endOrder<=length; endOrder+=10000) {
						try {
							pairHHI.appendArray(chrom, strand, currloaderPairHHI.getElement(chrom, strand, startOrder, endOrder));
						} catch (Exception e) {
							e.printStackTrace();
						}
						startOrder = endOrder;
					}
					// append the rest part of the hits into the hitcache
					if (currloaderPairHHI.getLength(chrom, strand)>startOrder)
						try {
							pairHHI.appendArray(chrom, strand, currloaderPairHHI.getElement(chrom, strand, startOrder, length));
						} catch (Exception e) {
							e.printStackTrace();
						}
				}
			
			// close the hitloader
			currLoader.close();
		}
		
		// Sort each hitInfo dataset
		readHHI.sortByReference();
		pairHHI.sortByReference();
	}
	
	//Accessors
	public double getHitCount(){return(totalHits);}
	public double getHitCountPos(){return(totalHitsPos);}
	public double getHitCountNeg(){return(totalHitsNeg);}
	public double getHitPositionCount(){return(uniqueHits);}
	public double getPairCount(){return(totalPairs);}
	public double getUniquePairCount(){return(uniquePairs);}
	public Genome getGenome(){return gen;}
	public void setGenome(Genome g){gen=g;}
	
	/**
	 * Set a new reference for the readHHI and sort the hits by this reference
	 * @param referenceString
	 */
	public void setReadReference(String referenceString) {
		readReference = referenceString;
		readHHI.setReference(referenceString);
		readHHI.sortByReference();
	}
	
	/**
	 * Set a new reference for the pairHHI and sort the hits by this reference
	 * @param referenceString
	 */
	public void setPairReference(String referenceString) {
		pairReference = referenceString;
		pairHHI.setReference(referenceString);
		pairHHI.sortByReference();
	}
	
	/**
	 * close the readhit and hitpair dataset
	 */
	public void close() {
		readHHI.closeDataset();
		pairHHI.closeDataset();
		readHHI.closeFile();
		pairHHI.closeFile();
	}
	
	/**
	 * Update the number of hits after loading all hitLoaders
	 */
	private void updateHits() {
		// Reset all hits number
		totalHits = 0;		
		totalHitsPos = 0;	
		totalHitsNeg = 0;	
		uniqueHits = 0;	
		totalPairs = 0;	
		uniquePairs = 0;	
		
		// Update number of hits associated with reads
		for(String chr: gen.getChromList())
			for(int strand=0; strand<2; strand++) {
				double[] weight = readHHI.getElement(chr, strand, "weight");
				for(int i=0; i<weight.length; i++) {
					totalHits += weight[i];
					uniqueHits++;
					if(strand==0)
						totalHitsPos++;
					else
						totalHitsNeg++;
				}
			}
		
		// Update number of hits associated with hitpairs
		for(String chr: gen.getChromList())
			for(int strand=0; strand<2; strand++) {
				double[] pairWeight = pairHHI.getElement(chr, strand, "pairWeight");
				for(int i=0; i<pairWeight.length; i++) {
					totalPairs += pairWeight[i];
					uniquePairs++;
				}
			}
	}
	
	/**
	 * Generate a hashmap containing the frequencies of each fragment size,
	 * @return
	 * @author Jianyu Yang
	 */
	public HashMap<Integer, Integer> getFragSizeFrequency(){
		HashMap<Integer, Integer> frequency = new HashMap<Integer, Integer>();
		if(loadPairs && hasPairs) {
			for(String chr: gen.getChromList())
				for(int strand=0; strand<2; strand++) {
					double[] r1Pos = pairHHI.getElement(chr, strand, "r1Pos");
					double[] r2Pos = pairHHI.getElement(chr, strand, "r2Pos");
					for(int i=0; i<r1Pos.length; i++) {
						int size = Math.abs((int)(r2Pos[i] - r1Pos[i]));
						if(frequency.containsKey(size)) {
							int oldValue = frequency.get(size);
							int newValue = oldValue + 1;
							frequency.put(size, newValue);
						} else {
							frequency.put(size, 1);
						}
					}
				}
		}
		return frequency;
	}
	
	/**
	 * Load all bases in a region according to the 5' end, regardless of the strand
	 * @param r
	 * @return
	 */
	public List<StrandedBaseCount> getBases(Region r) {
		List<StrandedBaseCount> bases = new ArrayList<StrandedBaseCount>();
		bases.addAll(getStrandedBases(r, '+'));
		bases.addAll(getStrandedBases(r, '-'));
		return bases;
	}
	
	/**
	 * Load all bases in a region according to the 5' end on specific strand
	 * @param r
	 * @param strand
	 * @return
	 */
	public synchronized List<StrandedBaseCount> getStrandedBases(Region r, char strand){
		if(readReference != "5'end")
			setReadReference("5'end");
		
		List<StrandedBaseCount> bases = new ArrayList<StrandedBaseCount>();
		
		String chr = r.getChrom();
		if(gen.getChromList().contains(chr)) {
			int j = (strand=='+') ? 0 : 1;
			int length = readHHI.getLength(chr, j);
			if(length > 0) {
				int start_ind = binarySearch(chr, j, r.getStart(), false);
				int end_ind = binarySearch(chr, j, r.getEnd(), false);
				
				if( start_ind < 0 ) {start_ind = -start_ind -1;}
				if( end_ind < 0)	{end_ind  = -end_ind - 1;}
				
				try {
					while(start_ind>0 && readHHI.getReference(chr, j, start_ind-1) >= r.getStart()) {
						start_ind--;
					}
					while(end_ind < length && readHHI.getReference(chr, j, end_ind) <= r.getEnd()) {
						end_ind++;
					}
					bases.addAll(readHHI.getBase(chr, j, start_ind, end_ind));
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
		}
		
		return bases;
	}
	
	/**
	 * Load all hit pairs in the given region according to the r1Pos, regardless of strand
	 * @param r
	 * @return
	 */
	public List<StrandedPair> getPairs(Region r) {
		List<StrandedPair> pairs = new ArrayList<StrandedPair>();
		pairs.addAll(getPairsOnStrand(r,'+'));
		pairs.addAll(getPairsOnStrand(r,'-'));
		return pairs;
	}
	
	/**
	 * Load all hit pairs in the given region according to the r1Pos on specific strand
	 * @param r
	 * @param strand
	 * @return
	 */
	public synchronized List<StrandedPair> getPairsOnStrand(Region r, char strand){
		if(pairReference != "r1Pos")
			setPairReference("r1Pos");
		
		List<StrandedPair> pairs = new ArrayList<StrandedPair>();
		
		if(loadPairs && hasPairs) {
			String chr = r.getChrom();
			if(gen.getChromList().contains(chr)) {
				int j = (strand == '+') ? 0 : 1;
				int length = pairHHI.getLength(chr, j);
				if(length>0) {
					int start_ind = binarySearch(chr, j, r.getStart(), true);
					int end_ind = binarySearch(chr, j, r.getEnd(), true);
					
					if( start_ind < 0 ) {start_ind = -start_ind -1;}
					if( end_ind < 0)	{end_ind  = -end_ind - 1;}
					
					try {
						
						while(start_ind>0 && pairHHI.getReference(chr, j, start_ind-1) >= r.getStart()) {
							start_ind--;
						}
						while(end_ind < length && pairHHI.getReference(chr, j, end_ind) <= r.getEnd()) {
							end_ind++;
						}
						pairs.addAll(pairHHI.getPair(chr, j, start_ind, end_ind));
					} catch (Exception e) {
						e.printStackTrace();
						System.exit(1);
					}
				}
			}
		}
		
		return pairs;
	}
	
	/**
	 * load all hit pairs in the given region according to the pairMid, regardless of the strand
	 * @param r
	 * @return
	 */
	public synchronized List<StrandedPair> getPairsByMid(Region r) {
		if(pairReference!="pairMid")
			setPairReference("pairMid");
		
		List<StrandedPair> pairs = new ArrayList<StrandedPair>();
		
		if(loadPairs && hasPairs) {
			String chr = r.getChrom();
			if(gen.getChromList().contains(chr)) 
				for(int j=0; j<2; j++) {
					int length = pairHHI.getLength(chr, j);
					if(length>0) {
						int start_ind = binarySearch(chr, j, r.getStart(), true);
						int end_ind = binarySearch(chr, j, r.getEnd(), true);
						
						if( start_ind < 0 ) {start_ind = -start_ind -1;}
						if( end_ind < 0)	{end_ind  = -end_ind - 1;}
						
						try {
							while(start_ind>0 && pairHHI.getReference(chr, j, start_ind-1) >= r.getStart()) {
								start_ind--;
							}
							while(end_ind < length && pairHHI.getReference(chr, j, end_ind) <= r.getEnd()) {
								end_ind++;
							}
							pairs.addAll(pairHHI.getPair(chr, j, start_ind, end_ind));
						} catch (Exception e) {
							e.printStackTrace();
							System.exit(1);
						}
					}
				}
		}
		
		return pairs;		
	}
	
	/**
	 * Same algorithm as Arrays.binarySearch, just to fit HierarchicalHitInfo
	 * @param chr
	 * @param strand
	 * @param index
	 * @param isPair
	 * @return
	 */
	private int binarySearch(String chr, int strand, int index, boolean isPair) {
		HierarchicalHitInfo hhInfo;
		if(isPair)
			hhInfo = pairHHI;
		else
			hhInfo = readHHI;
		
		int begin = 0;
		int last = hhInfo.getLength(chr, strand)-1;
		int mid = 0;
		
		while (begin <= last) {
			try {
				mid = (begin + last) >>> 1;
				double midIndex = hhInfo.getReference(chr, strand, mid);
				final int r = Double.compare(midIndex, index);
				if(r == 0) {
					return mid;
				} else if (r > 0) {
					last = mid - 1;
				} else {
					// This gets the insertion point right on the last loop
					begin = ++mid;
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		return -mid - 1;
	}
	
	public static void main(String[] args) {
		File genFile = new File(args[0]);
		File bam1File = new File(args[1]);
		File bam2File = new File(args[2]);
		
		Genome gen = new Genome("Genome", genFile, true);
		ExptConfig ec = new ExptConfig(gen, args);
		
		HDF5HitLoader hl1 = new HDF5HitLoader(gen, bam1File, false, false, false, true);
		HDF5HitLoader hl2 = new HDF5HitLoader(gen, bam2File, false, false, false, true);
		
		List<HDF5HitLoader> hList = new ArrayList<HDF5HitLoader>() {
			{
				add(hl1);
				add(hl2);
			}
		};
		
		HDF5HitCache hc = new HDF5HitCache(false, true, ec, hList);
		
		long start = System.currentTimeMillis();
		for(int i=0; i<10000; i++)
			hc.getPairsByMid(new Region(gen, "I", 10000, 10200));
		long end = System.currentTimeMillis();
		System.err.println((end - start) + "ms");
		hc.close();
		
	}
	
}
