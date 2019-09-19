package org.seqcode.deepseq.experiments;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.ExtReadHit;
import org.seqcode.deepseq.ReadHit;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.hitloaders.HDF5HitLoader;
import org.seqcode.deepseq.hitloaders.HitLoader;
import org.seqcode.deepseq.utils.HierarchicalHitInfo;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;


public class HDF5HitCache implements HitCacheInterface{
	
	private HierarchicalHitInfo readHHI;
	private HierarchicalHitInfo pairHHI;
	private Collection<HitLoader> loaders;	// source of hits
	private File localCacheDir;
	private Genome gen;
	private ExptConfig econfig;
	protected String name;
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
	
	public HDF5HitCache(ExptConfig ec, Collection<HitLoader> hloaders, String sampleName) {
		this.econfig = ec;
		this.gen = econfig.getGenome();
		this.loaders = hloaders;
		this.name = sampleName;
		
		this.loadReads = econfig.loadReads;
		this.loadPairs = econfig.loadPairs;
		
		this.localCacheDir = new File("");

		initialize();
	}
	
	/**
	 * load hit info of 		
		System.out.println(name);
		System.out.println(localCacheDir.getAbsolutePath());
		all hitloaders to a single dataset and sort them by reference element
	 */
	private void initialize() {
		
		// Create the hitInfo dataset for the HDF5HitCache itself
		readHHI = new HierarchicalHitInfo(gen, localCacheDir.getAbsolutePath()+ "/" + name.replace(':', '_') + ".read.h5", false);
		pairHHI = new HierarchicalHitInfo(gen, localCacheDir.getAbsolutePath()+ "/" + name.replace(':', '_') + ".pair.h5", true);
		
		// Initialize
		readHHI.initializeHDF5();
		pairHHI.initializeHDF5();
		
		// Append all hits into the dataset and extract the info into the HDF5HitCache
		boolean needSort = loaders.size()>1;
		for(HitLoader currLoader: loaders) {
			HDF5HitLoader currHDF5Loader = (HDF5HitLoader)currLoader;
			if(!currHDF5Loader.isCache()) {needSort = true;}
			currLoader.sourceAllHits();
			
			// readHHI
			if(loadReads) {
				HierarchicalHitInfo currloaderReadHHI = currHDF5Loader.getReadInfo();
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
			}
			
			// pairHHI
			if(loadPairs) {
				HierarchicalHitInfo currloaderPairHHI = currHDF5Loader.getHitPairInfo();
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
				
				// close the hitloader and delete the file
				currHDF5Loader.close();
				currHDF5Loader.cleanup();;
			}
		}
		
		// Sort each hitInfo dataset if needed
		if(needSort) {
			readHHI.sortByReference();
			pairHHI.sortByReference();
		}
		
		// Update hit number info
		updateHits();
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
	 * delete the readhit and hitpair dataset
	 */
	public void close() {
		readHHI.closeDataset();
		readHHI.closeFile();
		pairHHI.closeDataset();
		pairHHI.closeFile();
		if(!econfig.getKeepHDF5()) {
			readHHI.deleteFile();
			pairHHI.deleteFile();
		}
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
				if(readHHI.getLength(chr, strand)>0) {
					double[] weight = readHHI.getElement(chr, strand, "weight");
					for(int i=0; i<weight.length; i++) {
						totalHits += weight[i];
						uniqueHits++;
						if(strand==0)
							totalHitsPos++;
						else
							totalHitsNeg++;
					}
					weight = null;
				}
			}
		
		// Update number of hits associated with hitpairs
		for(String chr: gen.getChromList())
			for(int strand=0; strand<2; strand++) {
				if(pairHHI.getLength(chr, strand)>0) {
					double[] pairWeight = pairHHI.getElement(chr, strand, "pairWeight");
					for(int i=0; i<pairWeight.length; i++) {
						totalPairs += pairWeight[i];
						uniquePairs++;
					}
					pairWeight = null;
				}
			}
	}
	
	/**
	 * Simple count correction with a scaling factor and a floor of one. 
	 * Beware: only works if all reads are loaded.
	 * @param perBaseScaling float threshold
	 */
	public void linearCountCorrection(float perBaseScaling){
		if(perBaseScaling<1)
			System.err.println("linearCountCorrection: perBaseScaling is less than 1 - makes no sense to scale");
		else{
			for(String chr: gen.getChromList())
				for(int strand=0; strand<2; strand++)
					if(readHHI.getLength(chr, strand)>0) {
						double[] weights = readHHI.getElement(chr, strand, "weights");
						for( int i=0; i<weights.length; i++)
							if(weights[i] > 0) {
								weights[i] = weights[i] / perBaseScaling;
								if(weights[i] < 1)
									weights[i] = 1;
							}
						
						// write back the weights
						try {
							readHHI.updateElement(chr, strand, "weights", weights);
						} catch (Exception e) {
							e.printStackTrace();
						}
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
					if(pairHHI.getLength(chr, strand)>0) {
						double[] r1Pos = pairHHI.getElement(chr, strand, "r1Pos");
						double[] r2Pos = pairHHI.getElement(chr, strand, "r2Pos");
						for(int i=0; i<r1Pos.length; i++) {
							int size = Math.abs((int)(r2Pos[i] - r1Pos[i]));
							if(frequency.containsKey(size)) {
								frequency.put(size, frequency.get(size) + 1);
							} else {
								frequency.put(size, 1);
							}
						}
						r1Pos = null;
						r2Pos = null;
					}
				}
		}
		return frequency;
	}
	
	/**
	 * Sum of all hit weights in a region
	 * @param r Region
	 * @return count float
	 */
	public float countHits(Region r) {
		float count = 0;
		count += countStrandedBases(r, '+');
		count += countStrandedBases(r, '-');
		return count;
	}
	
	/**
	 * Sum of hit weights in one strand of a region
	 * @param r Region, strand char
	 * @return count float
	 */
	public synchronized float countStrandedBases(Region r, char strand) {
		// make sure the reference element of the read hit is 5' end
		if(readReference != "5'end")
			setReadReference("5'end");
		
		float count = 0;
		
		String chr = r.getChrom();
		if(gen.getChromList().contains(chr)) {
			int j = (strand=='+') ? 0 : 1;
			try {
				int[] ind = getIndexInRegion(chr, j, r, false);
				double[] weight = readHHI.getElement(chr, j, "5'end", ind[0], ind[1]);
				for(int i=0; i<weight.length; i++)
					count += weight[i];
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		return count;
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
			try {
				int[] ind = getIndexInRegion(chr, j, r, false);
				if(ind[1] > ind[0])
					bases.addAll(readHHI.getBase(chr, j, ind[0], ind[1]));
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
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
				try {
					int[] ind = getIndexInRegion(chr, j, r, true);
					if(ind[1] > ind[0])
						pairs.addAll(pairHHI.getPair(chr, j, ind[0], ind[1]));
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(1);
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
					try {
						int[] ind = getIndexInRegion(chr, j, r, true);
						if(ind[1] > ind[0])
							pairs.addAll(pairHHI.getPair(chr, j, ind[0], ind[1]));
					} catch (Exception e) {
						e.printStackTrace();
						System.exit(1);
					}
				}
		}
		
		return pairs;		
	}
	
	/**
	 * Exports readhits in a given region
	 */
	public List<ReadHit> exportReadHits(Region r, int readLen){
		List<ReadHit> reghits = new ArrayList<ReadHit>();
		for(StrandedBaseCount sbc: getBases(r)) {
    		int start = sbc.getStrand() == '+' ? sbc.getCoordinate() : sbc.getCoordinate()-readLen+1;
    		int end = sbc.getStrand() == '+' ? sbc.getCoordinate()+readLen-1 : sbc.getCoordinate();
    		reghits.add(new ReadHit(r.getChrom(),start,end,sbc.getStrand(),sbc.getCount()));
		}
		return reghits;
	}
	
    public List<ExtReadHit> exportExtReadHits(Region r, int readLen, int startShift, int fivePrimeExt, int threePrimeExt){
    	List<ReadHit> readHits = exportReadHits(r,readLen);
    	List<ExtReadHit> extReadHits = new ArrayList<ExtReadHit>();
    	for(ReadHit rh : readHits){
    		extReadHits.add(new ExtReadHit(gen, rh, startShift, fivePrimeExt, threePrimeExt));
    	}
    	return extReadHits;
    }
	
    /**
     * Export all single-end hits to ReadHit objects
     * @param readLen
     * @return
     */
    public List<ReadHit> exportReadHits(int readLen){
		List<ReadHit> allhits = new ArrayList<ReadHit>();
		
		for(String chr : gen.getChromList()){
			Region chrom = new Region(gen, chr, 1, gen.getChromLength(chr));
			for(StrandedBaseCount sbc : getBases(chrom)){
				int start = sbc.getStrand()=='+' ? sbc.getCoordinate() : sbc.getCoordinate()-readLen+1;
				int end = sbc.getStrand()=='+' ? sbc.getCoordinate()+readLen-1 : sbc.getCoordinate();
				allhits.add(new ReadHit(chr, start, end, sbc.getStrand(), sbc.getCount()));
			}
		}
		return(allhits);
	}
	
	/**
	 * Return the index of the hits with reference element located in the given region
	 * @param chr
	 * @param strand
	 * @param r
	 * @param isPair
	 * @return
	 */
	private int[] getIndexInRegion(String chr, int strand, Region r, boolean isPair) {
		int start_ind; int end_ind;
		HierarchicalHitInfo hhInfo = isPair ? pairHHI : readHHI;
		int length = hhInfo.getLength(chr, strand);
		if(length>0) {
			start_ind = binarySearch(chr, strand, r.getStart(), true);
			end_ind = binarySearch(chr, strand, r.getEnd(), true);
			
			if( start_ind < 0 ) {start_ind = -start_ind -1;}
			if( end_ind < 0)	{end_ind  = -end_ind - 1;}
			
			try {
				while(start_ind>0 && hhInfo.getReference(chr, strand, start_ind-1) >= r.getStart()) {
					start_ind--;
				}
				while(end_ind < length && hhInfo.getReference(chr, strand, end_ind) <= r.getEnd()) {
					end_ind++;
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		} else {
			return new int[] {0, 0};
		}
		
		return new int[] {start_ind, end_ind};
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
		
		HDF5HitLoader hl1 = new HDF5HitLoader(gen, bam1File, true, true, true, true, false);
		HDF5HitLoader hl2 = new HDF5HitLoader(gen, bam2File, true, true, true, true, false);
		
		List<HitLoader> hList = new ArrayList<HitLoader>() {
			{
				add(hl1);
				add(hl2);
			}
		};
		
		long start = System.currentTimeMillis();
		HDF5HitCache hc = new HDF5HitCache(ec, hList, "test");
		long end = System.currentTimeMillis();
		System.err.println((end - start) + "ms");
		
		List<StrandedPair> pairList = hc.getPairsByMid(new Region(gen, "I", 1000000000, 1000000100));
		for (StrandedPair sp: pairList) {
			System.out.println(sp);
		}
		List<StrandedBaseCount> baseList = hc.getBases(new Region(gen, "M", 5000, 5200));
		for (StrandedBaseCount sbc: baseList) {
			System.out.println(sbc);
		}
		List<StrandedPair> pairList2 = hc.getPairs(new Region(gen, "I", 3000, 3100));
		for (int i  = 0; i  < 5; i ++) {
			System.out.println(pairList2.get(i));
		}

		hc.close();
		
	}
	
}
