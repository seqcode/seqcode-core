package edu.psu.compbio.seqcode.projects.multigps.experiments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;
import edu.psu.compbio.seqcode.projects.multigps.framework.BackgroundCollection;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.PoissonBackgroundModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.*;

/**
 * Sample represents a single experimental sample whose reads are sourced from one or more HitLoaders. 
 * This class acts as a cache for all alignment hits associated with a particular experiment. 
 * This combines functionality from DeepSeq and ReadCache in the old setup.
 * 
 * Hit alignments are loaded into two primative type 3D arrays -- fivePrimePos and fivePrimeCounts.
 * Details from old ReadCache:
 * The fivePrimes field for each chrom/strand will be distinct. 
 * Multiple reads mapped to the same bp position will be stored as counts.
 * This class basically stores the hits coming from the corresponding files. <br>
 * We have made use of an unusual convention for reducing running time purposes. <br>
 * The hits are basically being represented by 3 main fields: <tt>fivePrimes, hitCounts</tt>
 * and <tt>hitIDs</tt>. <br>
 * Each of these fields are 3D arrays where:  <br>
 * - the first dimension corresponds to the chromosome that a hit belongs to (based on
 * the mapping from a chromosome as a <tt>String</tt> to an integer via the 
 * <tt>chrom2ID</tt> map).    <br>
 * - the second dimension corresponds to the strand. 0 for '+' (Watson), 1 for '-' (Crick). <br>
 * - the third dimension contains information for a hit (e.g. its fivePrimes or counts)
 * 
 * @author mahony
 *
 */
public class Sample {

	private int index;
	private Collection<HitLoader> loaders;
	private Config config;
	private Genome gen;
	protected String name;
	protected double totalHits; //totalHits is the sum of alignment weights
	private int numChroms;
	protected BackgroundCollection perBaseBack=new BackgroundCollection();
	protected float maxReadsPerBP=-1;
	
	/**
	 * Five prime ends of the read hits. <br>
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the coordinates of the hits
	 */
	private int[][][] fivePrimePos=null;
	/**
	 * Sum of read hit weights that corresponds to the 5' position
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the number of hits at corresponding start position 
	 */
	private float[][][] fivePrimeCounts=null;
	
	private HashMap<String, Integer> chrom2ID=new HashMap<String,Integer>();
	private HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
	
	/**
	 * Constructor
	 * @param g Genome (can be null to estimate from data)
	 * @param name String
	 */
	public Sample(int index, Config c, String name, float perBaseReadMax){
		this.index = index;
		this.name = name;
		config = c;
		gen=c.getGenome();
		totalHits=0;
		numChroms = 0;
		loaders = new ArrayList<HitLoader>();
		maxReadsPerBP= perBaseReadMax;
	}

	//Accessors
	public int getIndex(){return index;}
	public Genome getGenome(){return(gen);}
	public String getName(){return name;}
	public double getHitCount(){return(totalHits);}
	public void setGenome(Genome g){gen=g;}
	public BackgroundCollection getPerBaseBackground(){return perBaseBack;}

	/**
	 * Add a HitLoader to the set
	 * @param h HitLoader
	 */
	public void addHitLoader(HitLoader h){loaders.add(h);}
	
	/**
	 * Load hits from the loaders. Store everything in the primitive arrays.
	 * Hits from all sources are added to lists before converting to the primitive arrays 
	 * to facilitate collecting multiple sources in a single experiment. 
	 */
	public void loadHits(){
		//These lists are temporary stores while collecting reads from all sources
		HashMap<String, ArrayList<Integer>[]> posList = new HashMap<String, ArrayList<Integer>[]>();
		HashMap<String, ArrayList<Float>[]> countsList = new HashMap<String, ArrayList<Float>[]>();
		
		for(HitLoader currLoader : loaders){
			//Get the reads
			currLoader.sourceReads();
		
			//Add the reads to the temporary stores
			for(String chr: currLoader.getFivePrimePositions().keySet()){
				if(!posList.containsKey(chr)){
					ArrayList<Integer>[] currIArrayList = new ArrayList[2];
					currIArrayList[0]=new ArrayList<Integer>();
					currIArrayList[1]=new ArrayList<Integer>();
					posList.put(chr, currIArrayList);
				}
				posList.get(chr)[0].addAll(currLoader.getFivePrimePositions().get(chr)[0]);
				posList.get(chr)[1].addAll(currLoader.getFivePrimePositions().get(chr)[1]);
			}
			for(String chr: currLoader.getFivePrimeCounts().keySet()){
				if(!countsList.containsKey(chr)){
					ArrayList<Float>[] currFArrayList = new ArrayList[2];
					currFArrayList[0]=new ArrayList<Float>();
					currFArrayList[1]=new ArrayList<Float>();
					countsList.put(chr, currFArrayList);
				}
				countsList.get(chr)[0].addAll(currLoader.getFivePrimeCounts().get(chr)[0]);
				countsList.get(chr)[1].addAll(currLoader.getFivePrimeCounts().get(chr)[1]);
			}
			
			//Reset loader to free memory
			currLoader.resetLoader();
		}
		
		//Make the primitive arrays (null genome is estimated here if necessary)
		populateArrays(posList, countsList);
		
		//Free memory
		posList=null;
		countsList=null;
		System.gc();
		
		//Initialize a per-base background model
		initializeBackground();
		
		//Enforce a per-base read limit 
		//maxReadsPerBP = 0 : poisson/gauss
		//maxReadsPerBP = -1 : global poisson
		//maxReadsPerBP > 0 : fixed
		if(config.doPoissonGaussWinPerBaseFiltering() || maxReadsPerBP==0){ //global poisson/gauss model
			capPerBaseCountWithPoissonGaussianFilter(10e-3, 20);
		}else{
			if(maxReadsPerBP == -1)
				maxReadsPerBP = perBaseBack.getMaxThreshold('.');
			capPerBaseCount(maxReadsPerBP);
		}
		initializeBackground(); //Reinitialize given updated hit count (again - just the per-base background model)
	}
	
	
	/**
	 * Load all base counts in a region, regardless of strand
	 * @param r Region
	 * @return List of StrandedBaseCounts
	 */
	public List<StrandedBaseCount> getUnstrandedBases(Region r) {
		List<StrandedBaseCount> bases = new ArrayList<StrandedBaseCount>();
		bases.addAll(getStrandedBases(r,'+'));
		bases.addAll(getStrandedBases(r,'-'));
		return bases;
	}
	/**
	 * Loads hits in the region
	 * @param r Region
	 * @return List of StrandedBaseCounts
	 */
	public List<StrandedBaseCount> getStrandedBases(Region r, char strand) {
		List<StrandedBaseCount> bases = new ArrayList<StrandedBaseCount>();
		String chr = r.getChrom();
		int chrID = chrom2ID.get(chr);
		int j = (strand=='+') ? 0 : 1;
		if(fivePrimePos[chrID][j] != null){
			int[] tempStarts = fivePrimePos[chrID][j];		
			if(tempStarts.length != 0) {
				int start_ind = Arrays.binarySearch(tempStarts, r.getStart());
				int end_ind   = Arrays.binarySearch(tempStarts, r.getEnd());
				
				if( start_ind < 0 ) { start_ind = -start_ind - 1; }
				if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
				
	            while (start_ind > 0 && tempStarts[start_ind - 1] >= r.getStart() ) {
	                start_ind--;
	            }
	            while (end_ind < tempStarts.length && tempStarts[end_ind] <= r.getEnd()) {
	                end_ind++;
	            }
				for(int k = start_ind; k < end_ind; k++) {
					bases.add(new StrandedBaseCount(strand, tempStarts[k], fivePrimeCounts[chrID][j][k]));
				}	
			}
		}
		return bases;
	}//end of getStrandedBases method
	
	/**
	 * Sum of all hit weights in a region
	 * @param r Region
	 * @return float 
	 */
	public float countHits(Region r) {
		float count = 0;
		for (StrandedBaseCount b: getUnstrandedBases(r)){
			count += b.getCount();
		}
		return count;
	}
	/**
	 * Sum of hit weights in one strand of a region
	 * @param r Region
	 * @return float 
	 */
    public float countStrandedBases(Region r, char strand) {
		String chr = r.getChrom();
		int chrID = chrom2ID.get(chr);
		int j = (strand=='+') ? 0 : 1;
		float count = 0;
		if(fivePrimePos[chrID][j] != null){
			int[] tempStarts = fivePrimePos[chrID][j];		
	        if(tempStarts.length != 0) {
				int start_ind = Arrays.binarySearch(tempStarts, r.getStart());
				int end_ind   = Arrays.binarySearch(tempStarts, r.getEnd());
				if( start_ind < 0 ) { start_ind = -start_ind - 1; }
				if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
				
	            while (start_ind > 0 && tempStarts[start_ind - 1] >= r.getStart() ) {
	                start_ind--;
	            }
	            while (end_ind < tempStarts.length && tempStarts[end_ind] <= r.getEnd()) {
	                end_ind++;
	            }
				for(int k = start_ind; k < end_ind; k++) {
	                count += fivePrimeCounts[chrID][j][k];
	            }
	        }
		}
        return count;
    }

	/**
	 * Gets the stranded count of all hits (of all chromosomes) for the specified strand
	 * @param strand 
	 * @return
	 */
	protected double getStrandedTotalCount(char strand) {
		int strandInd = strand == '+' ? 0 : 1;
		double count = 0;
		for(int i = 0; i < fivePrimeCounts.length; i++) {
			if(fivePrimeCounts[i][strandInd]!=null){
				float[] hitCountsTemp = fivePrimeCounts[i][strandInd];
				for(float el:hitCountsTemp)
					count += (double)el;
			}
		}
		return count;
	}//end of getStrandedTotalCount method
	
	/**
	 * Converts lists of Integers to integer arrays, deletes the lists for saving memory
	 * Sorts array elements by position
	 */
	protected void populateArrays(HashMap<String, ArrayList<Integer>[]> posList, HashMap<String, ArrayList<Float>[]> countsList) {
		//Estimate the genome if none was provided
		if(gen==null){
			gen = estimateGenome(posList);
		}
		
		//Initialize chromosome name to id maps
		numChroms=0;
		for(String chr : gen.getChromList()){
			chrom2ID.put(chr, numChroms);
			id2Chrom.put(numChroms, chr);
			numChroms++;
		}
		
		//Initialize the data structures
		fivePrimePos  = new int[numChroms][2][];
		fivePrimeCounts = new float[numChroms][2][];
		
		//Copy over the data
		for(String chr : gen.getChromList()){
			if(posList.containsKey(chr)){
				for(int j = 0; j < posList.get(chr).length; j++)
					fivePrimePos[chrom2ID.get(chr)][j] = list2int(posList.get(chr)[j]);
			}else{
				fivePrimePos[chrom2ID.get(chr)][0]=null;
				fivePrimePos[chrom2ID.get(chr)][1]=null;
			}
		}
		for(String chr : gen.getChromList()){
			if(countsList.containsKey(chr)){
				for(int j = 0; j < countsList.get(chr).length; j++)
					fivePrimeCounts[chrom2ID.get(chr)][j] = list2float(countsList.get(chr)[j]);
			}else{
				fivePrimeCounts[chrom2ID.get(chr)][0]=null;
				fivePrimeCounts[chrom2ID.get(chr)][1]=null;
			}
		}
		//Sort the arrays 
		for(int i = 0; i < fivePrimePos.length; i++) {  // chr
			for(int j = 0; j < fivePrimePos[i].length; j++) { // strand
				if(fivePrimePos[i][j]!=null && fivePrimeCounts[i][j]!=null){
					int[] inds = StatUtil.findSort(fivePrimePos[i][j]);
					fivePrimeCounts[i][j] = StatUtil.permute(fivePrimeCounts[i][j], inds);
				}
			}
		}
		
		//Collapse duplicate positions
		//Testing if the finished set of positions contains duplicates
		for(int i = 0; i < fivePrimePos.length; i++){
			for(int j = 0; j < fivePrimePos[i].length; j++){
				if(fivePrimePos[i][j]!=null && fivePrimePos[i][j].length>0){
					int uniquePos=1;
					for(int k = 0; k < fivePrimePos[i][j].length-1; k++)
						if(fivePrimePos[i][j][k+1]!=fivePrimePos[i][j][k]){uniquePos++;}
							
					int[] tmpPos = new int[uniquePos];
					float[] tmpCnt = new float[uniquePos];
					for(int x=0; x<uniquePos; x++){tmpCnt[x]=0;}
					int x=0;
					tmpPos[x] = fivePrimePos[i][j][0];
					tmpCnt[x] += fivePrimeCounts[i][j][0];
					for(int k = 1; k < fivePrimePos[i][j].length; k++){
						if(fivePrimePos[i][j][k]!=fivePrimePos[i][j][k-1]){x++;}
						tmpPos[x] = fivePrimePos[i][j][k];
						tmpCnt[x] += fivePrimeCounts[i][j][k];
					}
					fivePrimeCounts[i][j] = tmpCnt;
					fivePrimePos[i][j] = tmpPos;
				}
			}
		}
		
		updateTotalHits();
		
		//Testing if the finished set of positions contains duplicates
		/*for(int i = 0; i < fivePrimePos.length; i++)
			for(int j = 0; j < fivePrimePos[i].length; j++)
				for(int k = 0; k < fivePrimePos[i][j].length-1; k++){
					if(fivePrimePos[i][j][k+1]==fivePrimePos[i][j][k])
						System.err.println("Duplicate position:\t"+i+"\t"+j+"\t"+k);
					else if(fivePrimeCounts[i][j][k]>1)
						System.err.println("Multi-position:\t"+i+"\t"+j+"\t"+fivePrimePos[i][j][k]+"\t"+fivePrimeCounts[i][j][k]);
				}
		*/
	}//end of populateArrays method
	
	/**
	 * Estimate a genome from the observed read positions that are collected into the list
	 * @param posList HashMap indexed by chr containing ArrayLists of hit positions
	 * @return Genome
	 */
	protected Genome estimateGenome(HashMap<String, ArrayList<Integer>[]> posList){
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		for(String c : posList.keySet()){
			int max = 0;
			for(int j=0; j<posList.get(c).length; j++){
				for(Integer i : posList.get(c)[j]){
					if(i>max)
						max=i;
				}
			}
			chrLenMap.put(c, max);
		}
		Genome g =new Genome("Genome", chrLenMap);
		return(g);
	}

	/**
	 * Recount hit weights
	 */
	protected void updateTotalHits(){
		totalHits = 0.0;
		for(int i = 0; i < fivePrimeCounts.length; i++)
			for(int j = 0; j < fivePrimeCounts[i].length; j++)
				if(fivePrimeCounts[i][j]!=null)
					for(int k = 0; k < fivePrimeCounts[i][j].length; k++)
						totalHits += fivePrimeCounts[i][j][k];
	}
	
	/**
	 * Initialize the per-base background model
	 *  
	 */
	protected void initializeBackground(){
		perBaseBack=new BackgroundCollection();
		perBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, config.getPerBaseLogConf(), getHitCount(), config.getGenome().getGenomeLength(), config.getMappableGenomeProp(), 1, '.', 1, true));
	}
	
	/**
	 * Simple convertor
	 * @param list
	 * @return
	 */
	protected int[] list2int(List<Integer> list) {
		int[] out = new int[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	/**
	 * Simple convertor
	 * @param list
	 * @return
	 */
	protected float[] list2float(List<Float> list) {
		float[] out = new float[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	
	/**
	 * Enforces a per-base weight threshold
	 * @param maxReadperBP float threshold
	 */
	public void capPerBaseCount(float maxReadperBP){
		for(int i = 0; i < fivePrimeCounts.length; i++)
			for(int j = 0; j < fivePrimeCounts[i].length; j++)
				if(fivePrimeCounts[i][j]!=null)
					for(int k = 0; k < fivePrimeCounts[i][j].length; k++)
						if (fivePrimeCounts[i][j][k] > maxReadperBP){
							//System.err.println("Capping "+fivePrimeCounts[i][j][k]+" to "+maxReadperBP);
							fivePrimeCounts[i][j][k] = maxReadperBP;
						}
						
		updateTotalHits();
	}
	
	/**
	 * Reset duplicate reads that pass Poisson threshold. 
	 * The Poisson lambda parameter is calculated by an Gaussian average
	 * that puts more weight for nearby bases (same chrom, same strand)
	 */
	public void capPerBaseCountWithPoissonGaussianFilter(double threshold, int width){
        double g[] = new double[width*4+1];
		NormalDistribution gaussianDist = new NormalDistribution(0, width*width);
		for (int i=0;i<g.length;i++)
			g[i]=gaussianDist.calcProbability((double)i);
		
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
			
		for(int i = 0; i < fivePrimeCounts.length; i++)
			for(int j = 0; j < fivePrimeCounts[i].length; j++){
				float counts[] = fivePrimeCounts[i][j];
				int pos[] = fivePrimePos[i][j]; 
				if(counts!=null){
					for(int k = 0; k < counts.length; k++){
						int posK = pos[k]; 
						double sum = 0;
						for (int x=1;x<=width*4;x++){		// at most extend out 250 idx
							if (k+x>=counts.length|| pos[k+x]-posK>width*4)
								break;
							sum += counts[k+x]*g[pos[k+x]-posK];
						}
						for (int x=1;x<=width*4;x++){		// at most extend out 250 idx
							if (k-x<0 || posK-pos[k-x]>width*4)
								break;
							sum += counts[k-x]*g[posK-pos[k-x]];
						}
						sum = sum/(1-g[0]);				// exclude this position for evaluation
						
						double countThres=0;
						P.setMean(sum);
						double pvalue=1;
						for(int b=1; pvalue>threshold; b++){
							pvalue=1-P.cdf(b);	//p-value as the tail of Poisson
							countThres=b;
						}
						if (counts[k] > Math.max(1,countThres))
							counts[k] = (float) Math.max(1,countThres);					
					}
				}
			}
		updateTotalHits();
	}

}
