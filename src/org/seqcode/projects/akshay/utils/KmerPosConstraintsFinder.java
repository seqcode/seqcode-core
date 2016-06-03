package org.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import java.util.Random;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.gse.utils.sequence.SequenceUtils;

/**
 * Given two sets of locations, the following code identifies kmer-kmer distances that are significantly found in the positive set over the negative set.
 * Usually used to analyse the k-mer features from an SVM analysis
 * 
 * Some hard coded constants in this code :-
 * 		Prints only those k-mer constraints whose prop of occurrence in the positive set is greater than 0.2
 * 											(AND)
 * 		Prints only those k-mer constraints whose prop of occurrence is greater in the pos set compared to neg set
 * 											(AND)
 * 		Prints only those k-mer constraints whose prop occurence zscore is greater than 5
 * 											(AND)
 * 		Prints only thise k-mer constraints whose p-value as calulate by smirnov test is less than or euwalt to 0.005
 * 
 * @author akshaykakumanu
 *
 */
public class KmerPosConstraintsFinder {
	
	public int k;          // Length of the k-mer
	public SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	public int winSize;    // The size of the window around the peak-pair location to consider for mapping the location into k-mer space
	public Genome gen;
	public List<Point> pos_points;      // List of points in the positive set
	public List<Region> pos_regions;    // positive points extended to regions with winSize as window size
	public List<Point> neg_points;     // List of points in the negative set
	public List<Region> neg_regions;   // negative points extended to regions with winSize as window size
	
	// Map matrices are 3-d int matrices with the following schema: -
	// 1st level :- kmer id
	// 2nd level :- region number in the pos/neg set
	// 3rd level :- position in the the 200bp window
	// value :- kmer id if ther is a match at that particular position/location. Else 0.
	// int[numk][neg_regions.size()][winSize]
	public int[][][] negMapMatrix;      
	public int[][][] posMapMatrix;
	
	// Pair matrices are 3-d int matrices with the following schema :-
	// 1st level :- region number in the pos/neg set
	// 2nd level :- kmer id (x)
	// 3rd level :- kmer id (y)
	// value :- min distance between kmer-x and kmer-y
	public int[][][] negPairMatrix; 
	public int[][][] posPairMatrix;
	
	// All the following have the following schema
	// 1st level :- kmer id (x)
	// 2nd level :- kmer id (y)
	// 3rd level :- distance beteween x and y
	// value :- look at the name
	public double[][][] pvalues;
	public double[][][] mean_pos;
	public double[][][] mean_neg;
	public double[][][] std_neg;
	public double[][][] z_scores;
	
	
	public static int sample_size = 100; // size of the sub sample
	public static int num_samples = 100; // number of sub somples
	public static double pvalue_c_level = 1.73;
	
	public KmerPosConstraintsFinder(Genome g, int win, int ksize) {
		this.gen = g;
		this.winSize = win;
		this.k = ksize;
	}
	
	
	// Settors
	
	/**
	 * loads points and regions for given filenames into pos and neg lists respectively
	 * @param posPeaksFileName
	 * @param negPeaksFileName
	 */
	public void setPointsAndRegions(String posPeaksFileName, String negPeaksFileName){
		this.pos_points = RegionFileUtilities.loadPeaksFromPeakFile(gen, posPeaksFileName, winSize);
		this.neg_points = RegionFileUtilities.loadPeaksFromPeakFile(gen, negPeaksFileName, winSize);
		this.pos_regions = RegionFileUtilities.loadRegionsFromPeakFile(gen, posPeaksFileName, winSize);
		this.neg_regions = RegionFileUtilities.loadRegionsFromPeakFile(gen, negPeaksFileName, winSize);
	}
	
	/**
	 * sets the map matrices both pos and neg
	 * @param SeqPathFile
	 */
	public void setMapMatrices(String SeqPathFile){
		int numk = (int)Math.pow(4, k);
		this.negMapMatrix = new int[numk][neg_regions.size()][winSize];
		this.posMapMatrix = new int[numk][pos_regions.size()][winSize];
		seqgen.useCache(true);
		seqgen.setGenomePath(SeqPathFile);
		
		for(int i=0; i<pos_regions.size(); i++){
			String seq = seqgen.execute(pos_regions.get(i)).toUpperCase();
			if(seq.contains("N"))
				continue;
			
			for(int j=0; j<(seq.length()-k+1); j++){
				String currK = seq.substring(j, j+k);
				String revCurrK = SequenceUtils.reverseComplement(currK);
				int currKInt = RegionFileUtilities.seq2int(currK);
				int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
				int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
				posMapMatrix[kmer][i][j] = kmer;
			}
		}
		
		for(int i=0; i<neg_regions.size(); i++){
			String seq = seqgen.execute(neg_regions.get(i)).toUpperCase();
			if(seq.contains("N"))
				continue;
			
			for(int j=0; j<(seq.length()-k+1); j++){
				String currK = seq.substring(j, j+k);
				String revCurrK = SequenceUtils.reverseComplement(currK);
				int currKInt = RegionFileUtilities.seq2int(currK);
				int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
				int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
				negMapMatrix[kmer][i][j] = kmer;
			}
		}
	}
	
	
	/**
	 * setting pair matrices by calling a getPairMatrix private methods (see that method for more details)
	 */
	public void setPairMatrices(){
		this.posPairMatrix = this.getPairMatrix(this.posMapMatrix);
		this.negPairMatrix = this.getPairMatrix(this.negMapMatrix);
	}
	
	/**
	 * Generates a list of unique integers from 0 to length, the size of the list the size of the sub sample
	 * @param length
	 * @return
	 */
	public List<Integer> getRandomIndex(int length){
		List<Integer> randIndexes = new ArrayList<Integer>();
		Random rand = new Random();
		int r = 0;
		do{
			do{
				r = rand.nextInt(length);
			}while(randIndexes.contains(r));
			randIndexes.add(r);
			
		}while(randIndexes.size() < KmerPosConstraintsFinder.sample_size);
		
		return randIndexes;
	}
	
	
	
	
	//Getters
	
	public int[][][] getPosMapMat(){return this.posMapMatrix;}
	public int[][][] getNegMapMat(){return this.negMapMatrix;}
	
	/**
	 * given a mapmat and a pairmat and other input params, this method returns the proportion 
	 * of the number of times a given distance off set is seen in the given list of ranInds
	 * @param mapmat
	 * @param pairmat
	 * @param kmerXind
	 * @param kmerYind
	 * @param randIndexes
	 * @return
	 */
	public double[] getSubMatricies(int[][][] mapmat, int[][][] pairmat, int kmerXind, int kmerYind, int[] randIndexes){
		double[] retPair = new double[this.winSize];
		for(int d=0; d<this.winSize; d++){
			int c = 0;
			for(int l=0; l<randIndexes.length; l++){
				if(pairmat[randIndexes[l]][kmerXind][kmerYind]  == d+1){
						c++;
				}
			}
			retPair[d]= (float)c/randIndexes.length;
		}
				
		return retPair;
	}
	
	/**
	 * Given a mapmat this generates a pairmat for that mapmat
	 * @param mapmat
	 * @return
	 */
	public int[][][] getPairMatrix(int[][][] mapmat){
		int numk = (int)Math.pow(4, k);
		int[][][] ret;
		ret =  new int[mapmat[0].length][numk][numk];
		for(int i=0; i<numk; i++){ // over all kmers (x)
			int irev = RegionFileUtilities.seq2int(SequenceUtils.reverseComplement(RegionFileUtilities.int2seq(i, k))); // Id of the reverse complement k-mer
			if(i>irev) // of the 2 complements, only consider the smaller one
				continue;
			for(int j=i; j<numk; j++){   // over all kmers (y), such that x<=y 
				int jrev = RegionFileUtilities.seq2int(SequenceUtils.reverseComplement(RegionFileUtilities.int2seq(j, k)));
				if(j>jrev) // of the 2 complements, only consider the smaller one
					continue;
				
				int kmerXind = i; // (x)
				int kmerYind =  j; // (y), such that x<=y and x and y are the chosen complements (smaller ones)
				for(int l=0; l< mapmat[0].length; l++){ // over all locations 
					int[] Xmap = mapmat[kmerXind][l];     // Map of (x) for location l
					int[] Ymap = mapmat[kmerYind][l];     // Map of (y) for location l
					int[] XYmap = new int[Xmap.length];   // Adding both maps ..
					for(int k=0; k<XYmap.length; k++){
						XYmap[k] = Xmap[k]+Ymap[k];
					}
					
					// the following code is to calculate the min distance between (x) and (y) in the location l
					
					int distance = 0;   // this will eventually store the min distance                       
					int currMinDistance = Integer.MAX_VALUE; // the value of the current min distance
					boolean counting = false;  // are you counting ?? No
					int leftInd = 0; // value of the leftIND, set to 0 to start with (stores the kmer id of the left ind)
					
					for(int k=0; k<XYmap.length; k++){ // for each locaiton in the 200bp window 
						if(!counting && XYmap[k] != 0){ // if you are not counting and (x) or (y) match to the current position 
							counting = true; // start counting now
							leftInd = XYmap[k]; // store value of the kmer id
							distance = 1; // set distance to 1
						}
						else if(counting && XYmap[k] == 0){ // if counting and no match, continue counting
							distance++;
						}
						else if(counting && XYmap[k] !=0 && XYmap[k] == leftInd){ // if counting and match and match is same as left ind
							if(kmerXind == kmerYind){ // if (x) and (y) are same
								if(distance < currMinDistance){ // if the calc distance less than curr min distance
									currMinDistance = distance; 
								}
								distance = 1; // reset distance to 1
								leftInd = XYmap[k]; // store left ind again
							}
							else{ // if (x) and (y) don't match reset counting without storing anything
								distance =1; // reset distane 
								leftInd = XYmap[k]; // store x ind
							}
						}
						else if(counting && XYmap[k] !=0 && XYmap[k] != leftInd){ // if counting and match and different from left ind match 
							if(distance < currMinDistance){ // compare with curr min, if less overwrite
								currMinDistance = distance;
							}
							distance = 1; // reset
							leftInd = XYmap[k]; // store left ind
						}
					}
					
					
					if(currMinDistance == Integer.MAX_VALUE){ // if no min distance found, store it as 200.. will remove later 
						ret[l][kmerXind][kmerYind] = 200;
					}else{
						ret[l][kmerXind][kmerYind] = currMinDistance ;
					}
					
				}
				
				//for(int l=0; l<this.winSize; l++){
				//	ret[l][kmerXind][kmerYind] = ret[l][kmerXind][kmerYind]/mapmat[0].length;
				//}
			
			}
		}
		
		return ret;
		
	}
	
	
	
	
	/**
	 * 
	 * @param kmerSet
	 */
	public void setPvalues(List<String> kmerSet){
		int numk = (int)Math.pow(4, k);
		this.pvalues = new double[numk][numk][this.winSize]; 
		this.mean_neg = new double[numk][numk][this.winSize]; 
		this.mean_pos = new double[numk][numk][this.winSize];
		this.std_neg = new double[numk][numk][this.winSize];
		this.z_scores = new double[numk][numk][this.winSize];
		
		//double[][][][] PosPDF = new double[this.winSize][numk][numk][KmerPosConstraintsFinder.num_samples];
		//double[][][][] NegPDF = new double[this.winSize][numk][numk][KmerPosConstraintsFinder.num_samples];
		
		// 1st level :- sub sample id
		// 2nd level :- sample number in the given subsample
		// value :- index in the original data
		int[][] PosRandIndexes = new int[KmerPosConstraintsFinder.num_samples][KmerPosConstraintsFinder.sample_size]; 
		int[][] NegRandIndexes = new int[KmerPosConstraintsFinder.num_samples][KmerPosConstraintsFinder.sample_size];
		for(int itr = 0; itr< KmerPosConstraintsFinder.num_samples; itr++){ // over the number of subsamples
			List<Integer> tempP = this.getRandomIndex(this.pos_regions.size()); // get random index list 
			for(int r=0; r< tempP.size(); r++){ // over the no of samples in a given subsample
				PosRandIndexes[itr][r] = tempP.get(r);
			}
			
			List<Integer> tempN = this.getRandomIndex(this.neg_points.size());
			for(int r=0; r<tempN.size(); r++){
				NegRandIndexes[itr][r] = tempN.get(r);
			}
		}
		
		// 1st level :- distance off set
		// 2nd level :- subsample number
		// value :- proportion of sites having the offset
		double[][] PosPDF = new double[this.winSize][KmerPosConstraintsFinder.num_samples];
		double[][] NegPDF = new double[this.winSize][KmerPosConstraintsFinder.num_samples];
		
		for(int i=0; i<kmerSet.size(); i++){ // over all the given list of kmers (x)
			for(int j=i; j<kmerSet.size(); j++){ // over all the given list, such that (y) >= (x)
				int iInd = RegionFileUtilities.seq2int(kmerSet.get(i)); // geting the index of (x)
				int iRevInd = RegionFileUtilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(i))); // getting rev index of (x)
				int Xind = iInd < iRevInd ? iInd: iRevInd; // set the smaller one to (x)
				int jInd = RegionFileUtilities.seq2int(kmerSet.get(j));
				int jRevInd = RegionFileUtilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(j)));
				int Yind = jInd < jRevInd ? jInd : jRevInd;
				
				int kmerXind = Xind < Yind ? Xind: Yind; // redefine (x) again
				int kmerYind = Xind < Yind ? Yind : Xind; // redefine (y) again
				
				for( int itr=0; itr< KmerPosConstraintsFinder.num_samples; itr++){ // over number of subsamples
					// get the proportions for each subsample
					double[] tempPosPair = this.getSubMatricies(posMapMatrix, posPairMatrix, kmerXind, kmerYind, PosRandIndexes[itr]);
					double[] tempNegPair = this.getSubMatricies(negMapMatrix, negPairMatrix, kmerXind, kmerYind, NegRandIndexes[itr]);
					for(int d=0; d<this.winSize; d++){ // over all off sets
						PosPDF[d][itr] = tempPosPair[d]; // fill
						NegPDF[d][itr] = tempNegPair[d];
					}
				}
				
				for(int d=0; d<this.winSize; d++){ // over all offsets
					double maxD = computeD(PosPDF[d], NegPDF[d]); // compute maxD (see kolmogrov-smirnov test) for each offset
					double test = KmerPosConstraintsFinder.pvalue_c_level * (Math.sqrt(2/(float)KmerPosConstraintsFinder.num_samples));
					if(maxD > test){
						this.pvalues[kmerXind][kmerYind][d] = 0.005;
						this.mean_pos[kmerXind][kmerYind][d] = getMean(PosPDF[d]);
						this.mean_neg[kmerXind][kmerYind][d] = getMean(NegPDF[d]);
						this.std_neg[kmerXind][kmerYind][d] = this.getSTD(NegPDF[d]);
						this.z_scores[kmerXind][kmerYind][d] = (this.mean_pos[kmerXind][kmerYind][d] - this.mean_neg[kmerXind][kmerYind][d])/this.std_neg[kmerXind][kmerYind][d];
						
					}else{
						this.pvalues[kmerXind][kmerYind][d] = 1.0;
					}
				}
				
				
			}
		}
	}
	
	
	// Calculators and printers
	
	public void printSigFIConstrains(List<String> kmerSet){
		for(int i=0; i<kmerSet.size(); i++){
			for(int j=i; j<kmerSet.size(); j++){
				int iInd = RegionFileUtilities.seq2int(kmerSet.get(i));
				int iRevInd = RegionFileUtilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(i)));
				int Xind = iInd < iRevInd ? iInd: iRevInd;
				
				int jInd = RegionFileUtilities.seq2int(kmerSet.get(j));
				int jRevInd = RegionFileUtilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(j)));
				int Yind = jInd < jRevInd ? jInd : jRevInd;
				
				int kmerXind = Xind < Yind ? Xind: Yind;
				int kmerYind = Xind < Yind ? Yind : Xind;
				List<Integer> posCon = new ArrayList<Integer>(); 
				for(int d=0; d< this.winSize; d++){
					if(this.pvalues[kmerXind][kmerYind][d] == 0.005 && this.mean_pos[kmerXind][kmerYind][d] > this.mean_neg[kmerXind][kmerYind][d] && this.mean_pos[kmerXind][kmerYind][d] >0.2  && this.z_scores[kmerXind][kmerYind][d] >5.0 && d!= 199){
						posCon.add(d);
					}
				}
				
				if(posCon.size() > 0){
					for(int d : posCon){
						String zscore = Double.toString((this.mean_pos[kmerXind][kmerYind][d] - this.mean_neg[kmerXind][kmerYind][d])/this.std_neg[kmerXind][kmerYind][d]);
						String outPropPos = Double.toString(this.mean_pos[kmerXind][kmerYind][d]);
						String outPropNeg = Double.toString(this.mean_neg[kmerXind][kmerYind][d]);
						System.out.println(RegionFileUtilities.int2seq(kmerXind, k)+"-"+RegionFileUtilities.int2seq(kmerYind, k)+"\t"+Integer.toString(d+1)+"\t"+outPropPos+"\t"+outPropNeg+"\t"+zscore);
						
					}
					
				}
				
				
				
				
				
				//if(posCon.size() > 0){
				//	String outPOS = "";
				//	String outPropPos = "";
				//	String outPropNeg = "";
				//	String outZscore = "";
				//	for(int d :  posCon){
				//		outPOS = outPOS+Integer.toString(d)+":";
				//		double zscore = (this.mean_pos[kmerXind][kmerYind][d] - this.mean_neg[kmerXind][kmerYind][d])/this.std_neg[kmerXind][kmerYind][d];
				//		outPropPos = outPropPos + Double.toString(this.mean_pos[kmerXind][kmerYind][d])+":";
				//		outPropNeg = outPropNeg + Double.toString(this.mean_neg[kmerXind][kmerYind][d])+":";
				//		outZscore = outZscore + Double.toString(zscore)+":";
				//	}
					
				//	System.out.println(RegionFileUtilities.int2seq(kmerXind, k)+"-"+RegionFileUtilities.int2seq(kmerYind, k)+"\t"+outPOS.substring(0,outPOS.length()-1 )+"\t"+outPropPos.substring(0, outPropPos.length()-1)
				//			+"\t"+outPropNeg.substring(0, outPropNeg.length()-1)+"\t"+outZscore.substring(0, outZscore.length()-1));
				//}
				
				
			}
		}
	}
	

	public double computeD(double[] vecx, double[] vecy){
		double maxDis;
		double[] vecCumX = new double[101];
		double[] vecCumY = new double[101];
		for(int i=0; i<=100; i++){
			double val = 0+i*0.01;
			for(int j=0; j<vecx.length; j++){
				if(vecx[j] <= val){
					vecCumX[i]++;
				}
				if(vecy[j] <= val){
					vecCumY[i]++;
				}
			}
		}
		double[] distances = new double[100];
		for(int i=0; i<distances.length; i++){
			double dis = (vecCumX[i]/vecx.length) - (vecCumY[i]/vecy.length);
			distances[i] = dis>0? dis : -1*dis;
		}
		
		maxDis = distances[getMaxIndex(distances)];
		return maxDis;
	}
	
	public double getSTD(double[] vec){
		double ret = 0.0;
		double mean = getMean(vec);
		
		for(int i=0; i<vec.length; i++){
			ret = ret + (vec[i] - mean)*(vec[i] - mean);
		}
		
		ret = ret/vec.length;
		
		ret = Math.sqrt(ret);
		
		return ret;
	}
	
	public double getMean(double[] vec){
		double ret=0.0;
		
		for(int i=0; i<vec.length; i++){
			ret = ret + vec[i];
		}
		
		return ret/vec.length;
	}
	public int getMinIndex(double[] vec){
		int ret=0;
		double currMin = Double.MAX_VALUE;
		for(int i=0; i< vec.length; i++){
			if(vec[i] < currMin){
				currMin = vec[i];
				ret = i;
			}
		}
		return ret;
	}
	
	public int getMaxIndex(double[] vec){
		int ret=0;
		double currMax = Double.MIN_VALUE;
		for(int i=0; i<vec.length; i++){
			if(vec[i]> currMax){
				currMax = vec[i];
				ret = i;
			}
		}
		return ret;
	}
	
	
	// main methods follows 
	
	public static void main(String[] args) throws IOException, NotFoundException{
		ArgParser ap = new ArgParser(args);
		System.err.println("Usage:\n"
		);
		//reading selected kmers
		String kmerFile = ap.getKeyValue("kmers");
		File kmerF = new File(kmerFile);
		BufferedReader reader = new BufferedReader(new FileReader(kmerF));
		String line;
		List<String> kmers = new ArrayList<String>();
		while((line = reader.readLine()) != null){
			line = line.trim();
			line.replace("\n", "");
			kmers.add(line);
		}
		reader.close();
		
		String SeqFilePath = ap.getKeyValue("seq");
		int winSize = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue() : 200;
		Genome gen = null;
		if(ap.hasKey("species")){
			Pair<Species, Genome> pair = Args.parseGenome(args);
			if(pair != null){
				
				gen = pair.cdr();
			}
		}else{
			if(ap.hasKey("geninfo") || ap.hasKey("g")){
				//Make fake genome... chr lengths provided
				String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
				gen = new Genome("Genome", new File(fName), true);
			}else{
				gen = null;
			}
		}
		
		
		int ksize = ap.hasKey("ksize") ? new Integer(ap.getKeyValue("ksize")).intValue() : 4;
		
		String PosPeaksFile = ap.getKeyValue("peaksP");
		String NegPeaksFile = ap.getKeyValue("peaksN");
		
		double zscoreCutoff = Args.parseDouble(args, "zCut", 5.5);
		double fractionCutoff = Args.parseDouble(args, "frac", 2.0);
		
		KmerPosConstraintsFinder analysis = new KmerPosConstraintsFinder(gen, winSize, ksize);
		
		analysis.setPointsAndRegions(PosPeaksFile, NegPeaksFile);
		analysis.setMapMatrices(SeqFilePath);
		analysis.setPairMatrices();
		
		analysis.setPvalues(kmers);
		
		analysis.printSigFIConstrains(kmers);
		
		
	}
	
	
	

}
