package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;
import java.util.Random;

public class KmerPosConstraintsFinder {
	
	public int k;
	public SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	public int winSize;
	public Genome gen;
	public List<Point> pos_points;
	public List<Region> pos_regions;
	public List<Point> neg_points;
	public List<Region> neg_regions;
	public int[][][] negMapMatrix;
	public int[][][] posMapMatrix;
	public int[][][] negPairMatrix; // location, xkmer-id, ykmer-id, min distance
	public int[][][] posPairMatrix;
	
	public double[][][] pvalues;
	public double[][][] mean_pos;
	public double[][][] mean_neg;
	public double[][][] std_neg;
	
	public static int sample_size = 100;
	public static int num_samples = 100;
	public static double pvalue_c_level = 1.73;
	
	public KmerPosConstraintsFinder(Genome g, int win, int ksize) {
		this.gen = g;
		this.winSize = win;
		this.k = ksize;
	}
	
	
	// Settors
	public void setPointsAndRegions(String posPeaksFileName, String negPeaksFileName){
		this.pos_points = Utilities.loadPeaksFromPeakFile(gen, posPeaksFileName, winSize);
		this.neg_points = Utilities.loadPeaksFromPeakFile(gen, negPeaksFileName, winSize);
		this.pos_regions = Utilities.loadRegionsFromPeakFile(gen, posPeaksFileName, winSize);
		this.neg_regions = Utilities.loadRegionsFromPeakFile(gen, negPeaksFileName, winSize);
	}
	
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
				int currKInt = Utilities.seq2int(currK);
				int revCurrKInt = Utilities.seq2int(revCurrK);
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
				int currKInt = Utilities.seq2int(currK);
				int revCurrKInt = Utilities.seq2int(revCurrK);
				int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
				negMapMatrix[kmer][i][j] = kmer;
			}
		}
	}
	
	
	
	public void setPairMatrices(){
		this.posPairMatrix = this.getPairMatrix(this.posMapMatrix);
		this.negPairMatrix = this.getPairMatrix(this.negMapMatrix);
	}
	
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
	public double[] getSubMatricies(int[][][] mapmat, int[][][] pairmat, int kmerXind, int kmerYind, int[] randIndexes){
		double[] retPair = new double[this.winSize];
		//List<Integer> randIndexes = new ArrayList<Integer>();
		//Random rand = new Random();
		//int r = 0;
		//do{
		//	do{
		//		r = rand.nextInt(mapmat[0].length);
		//	}while(randIndexes.contains(r));
		//	randIndexes.add(r);
			
		//}while(randIndexes.size() < KmerPosConstraintsFinder.sample_size);
		
		
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
	
	public int[][][] getPairMatrix(int[][][] mapmat){
		int numk = (int)Math.pow(4, k);
		int[][][] ret;
		ret =  new int[mapmat[0].length][numk][numk];
		for(int i=0; i<numk; i++){
			int irev = Utilities.seq2int(SequenceUtils.reverseComplement(Utilities.int2seq(i, k)));
			if(i>irev)
				continue;
			for(int j=i; j<numk; j++){
				int jrev = Utilities.seq2int(SequenceUtils.reverseComplement(Utilities.int2seq(j, k)));
				if(j>jrev)
					continue;
				
				int kmerXind = i;
				int kmerYind =  j;
				for(int l=0; l< mapmat[0].length; l++){ // over all locations 
					int[] Xmap = mapmat[kmerXind][l];
					int[] Ymap = mapmat[kmerYind][l];
					int[] XYmap = new int[Xmap.length];
					for(int k=0; k<XYmap.length; k++){
						XYmap[k] = Xmap[k]+Ymap[k];
					}
					int distance = 0;
					int currMinDistance = Integer.MAX_VALUE;
					boolean counting = false;
					int leftInd = 0;
					for(int k=0; k<XYmap.length; k++){
						if(!counting && XYmap[k] != 0){
							counting = true;
							leftInd = XYmap[k];
							distance = 1;
						}
						else if(counting && XYmap[k] == 0){
							distance++;
						}
						else if(counting && XYmap[k] !=0 && XYmap[k] == leftInd){
							if(kmerXind == kmerYind){
								if(distance < currMinDistance){
									currMinDistance = distance;
								}
								distance = 1;
								leftInd = XYmap[k];
							}
							else{
								distance =1;
								leftInd = XYmap[k];
							}
						}
						else if(counting && XYmap[k] !=0 && XYmap[k] != leftInd){
							if(distance < currMinDistance){
								currMinDistance = distance;
							}
							distance = 1;
							leftInd = XYmap[k];
						}
					}
					if(currMinDistance == Integer.MAX_VALUE){
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
	
	
	
	

	public void setPvalues(List<String> kmerSet){
		int numk = (int)Math.pow(4, k);
		this.pvalues = new double[numk][numk][this.winSize];
		this.mean_neg = new double[numk][numk][this.winSize];
		this.mean_pos = new double[numk][numk][this.winSize];
		this.std_neg = new double[numk][numk][this.winSize];
		//double[][][][] PosPDF = new double[this.winSize][numk][numk][KmerPosConstraintsFinder.num_samples];
		//double[][][][] NegPDF = new double[this.winSize][numk][numk][KmerPosConstraintsFinder.num_samples];
		
		
		int[][] PosRandIndexes = new int[KmerPosConstraintsFinder.num_samples][KmerPosConstraintsFinder.sample_size];
		int[][] NegRandIndexes = new int[KmerPosConstraintsFinder.num_samples][KmerPosConstraintsFinder.sample_size];
		for(int itr = 0; itr< KmerPosConstraintsFinder.num_samples; itr++){
			List<Integer> tempP = this.getRandomIndex(this.pos_regions.size());
			for(int r=0; r< tempP.size(); r++){
				PosRandIndexes[itr][r] = tempP.get(r);
			}
			
			List<Integer> tempN = this.getRandomIndex(this.neg_points.size());
			for(int r=0; r<tempN.size(); r++){
				NegRandIndexes[itr][r] = tempN.get(r);
			}
		}
		
		double[][] PosPDF = new double[this.winSize][KmerPosConstraintsFinder.num_samples];
		double[][] NegPDF = new double[this.winSize][KmerPosConstraintsFinder.num_samples];
		
		for(int i=0; i<kmerSet.size(); i++){
			for(int j=i; j<kmerSet.size(); j++){
				int iInd = Utilities.seq2int(kmerSet.get(i));
				int iRevInd = Utilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(i)));
				int Xind = iInd < iRevInd ? iInd: iRevInd;
				
				int jInd = Utilities.seq2int(kmerSet.get(j));
				int jRevInd = Utilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(j)));
				int Yind = jInd < jRevInd ? jInd : jRevInd;
				
				int kmerXind = Xind < Yind ? Xind: Yind;
				int kmerYind = Xind < Yind ? Yind : Xind;
				
				for( int itr=0; itr< KmerPosConstraintsFinder.num_samples; itr++){
					double[] tempPosPair = this.getSubMatricies(posMapMatrix, posPairMatrix, kmerXind, kmerYind, PosRandIndexes[itr]);
					double[] tempNegPair = this.getSubMatricies(negMapMatrix, negPairMatrix, kmerXind, kmerYind, NegRandIndexes[itr]);
					for(int d=0; d<this.winSize; d++){
						PosPDF[d][itr] = tempPosPair[d];
						NegPDF[d][itr] = tempNegPair[d];
					}
				}
				
				for(int d=0; d<this.winSize; d++){
					double maxD = computeD(PosPDF[d], NegPDF[d]);
					double test = KmerPosConstraintsFinder.pvalue_c_level * (Math.sqrt(2/(float)KmerPosConstraintsFinder.num_samples));
					if(maxD > test){
						this.pvalues[kmerXind][kmerYind][d] = 0.005;
						this.mean_pos[kmerXind][kmerYind][d] = getMean(PosPDF[d]);
						this.mean_neg[kmerXind][kmerYind][d] = getMean(NegPDF[d]);
						this.std_neg[kmerXind][kmerYind][d] = this.getSTD(NegPDF[d]);
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
				int iInd = Utilities.seq2int(kmerSet.get(i));
				int iRevInd = Utilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(i)));
				int Xind = iInd < iRevInd ? iInd: iRevInd;
				
				int jInd = Utilities.seq2int(kmerSet.get(j));
				int jRevInd = Utilities.seq2int(SequenceUtils.reverseComplement(kmerSet.get(j)));
				int Yind = jInd < jRevInd ? jInd : jRevInd;
				
				int kmerXind = Xind < Yind ? Xind: Yind;
				int kmerYind = Xind < Yind ? Yind : Xind;
				List<Integer> posCon = new ArrayList<Integer>(); 
				for(int d=0; d< this.winSize; d++){
					if(this.pvalues[kmerXind][kmerYind][d] == 0.005 && this.mean_pos[kmerXind][kmerYind][d] > this.mean_neg[kmerXind][kmerYind][d] && this.mean_pos[kmerXind][kmerYind][d] > 0.1 && d!= 199){
						posCon.add(d);
					}
				}
				
				if(posCon.size() > 0){
					String outPOS = "";
					String outProp = "";
					String outZscore = "";
					for(int d :  posCon){
						outPOS = outPOS+Integer.toString(d)+":";
						double zscore = (this.mean_pos[kmerXind][kmerYind][d] - this.mean_neg[kmerXind][kmerYind][d])/this.std_neg[kmerXind][kmerYind][d];
						outProp = outProp + Double.toString(this.mean_pos[kmerXind][kmerYind][d])+":";
						outZscore = outZscore + Double.toString(zscore)+":";
					}
					
					System.out.println(Utilities.int2seq(kmerXind, k)+" - "+Utilities.int2seq(kmerYind, k)+"\t"+outPOS+"\t"+outProp+"\t"+outZscore);
				}
				
				
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
			Pair<Organism, Genome> pair = Args.parseGenome(args);
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
		
		KmerPosConstraintsFinder analysis = new KmerPosConstraintsFinder(gen, winSize, ksize);
		
		analysis.setPointsAndRegions(PosPeaksFile, NegPeaksFile);
		analysis.setMapMatrices(SeqFilePath);
		analysis.setPairMatrices();
		
		analysis.setPvalues(kmers);
		
		analysis.printSigFIConstrains(kmers);
		
		
	}
	
	
	

}
