package edu.psu.compbio.seqcode.projects.naomi.xoqualitycontrol;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig; 
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.projects.naomi.FeatureCountsLoader;

/**
 * ReadDistributionProfiler : calculates weighted standard deviation around the reference point as in Rohit thesis
 * 
 * input : reference point
 * 
 * @author naomi yamada
 */

public class ReadDistributionProfiler {
	protected FeatureCountsLoader featureCountsLoader;
		
	protected List<StrandedPoint> strandedPoints;
	protected List<StrandedRegion> strandedRegions;
	
	protected int fivePrimeShift = 6;
	protected int window = 1000;
	protected int iterations = 1000;	
	protected boolean printSampleComposite = false;
	
	protected Map<Sample,double[]> sampleComposite = new HashMap<Sample,double[]>();
	protected Map<Sample,double[]> sampleStandardDeviation = new HashMap<Sample,double[]>();
	protected Map<Sample,double[]> sampleStatsAndNull = new HashMap<Sample,double[]>();
	
	public ReadDistributionProfiler(FeatureCountsLoader fcloader){	
		featureCountsLoader = fcloader;
	}
	
	// setters
	public void turnOnPrintComposite(){printSampleComposite = true;}
		
	public void StandardDeviationFromExpAndNull(){
		
		sampleComposite = featureCountsLoader.sampleComposite();
		

		for (Sample sample : sampleComposite.keySet()){
			// calculate standard deviation from sample
			double[] composite = sampleComposite.get(sample);			
			sampleStandardDeviation.put(sample, computeWeightedStandardDeviation(composite));
			if (printSampleComposite==true){
				System.out.println(sample.getName());
				printArray(composite);
			}
		
			//shuffle the data points and calculate standard deviation from null	
			double [] arrNullSD = new double [iterations]; 
			for (int itr = 0 ; itr < iterations; itr++){				
				double[] shuffledComposite =  shuffleComposite(composite);		
				arrNullSD[itr] = computeWeightedStandardDeviation(shuffledComposite)[1];
			}
			
			//for test purpose
			System.out.println("sd from null distribution");
			printArray(arrNullSD);		
			
			double mu = TotalCounts(arrNullSD)/arrNullSD.length;
			double sd = computeStandardDeviation(arrNullSD);		
			double x = sampleStandardDeviation.get(sample)[1];
			
			double[] x2 = new double [1];
			x2[0] = x;
			
			//Mann-Whitney test
//			MannWhitneyUTest mw = new MannWhitneyUTest();
//			double mannWP = mw.mannWhitneyUTest(x2,arrNullSD);
//			System.out.println("MannWhitney p value is "+mannWP);
			
			double z_score = (x - mu)/sd;
			double p_val = Erf.erfc(Math.abs(z_score)/Math.sqrt(2));
			
			double [] statistics = new double [5];
			statistics[0] = p_val;
			if ((x < mu) && p_val <0.05){statistics[1]=1;}else{statistics[1] = 0;}
			statistics[2] = z_score;
			statistics[3] = mu;
			statistics[4] = sd;			
			
			sampleStatsAndNull.put(sample, statistics);
		}
	}
	
	// This is for test purpose
	public void StandardDeviationFromRandomReads(){
		
		for (Sample sample : sampleComposite.keySet()){
			double[] composite = sampleComposite.get(sample);	
			double [] randomReadsSD = new double[iterations];
			for (int itr = 0 ; itr <iterations; itr++){
				double[] randomReads = randomelyAssignReads(composite);
				randomReadsSD[itr] = computeWeightedStandardDeviation(randomReads)[1];
			}
			System.out.println("standard deviation from random reads");
			printArray(randomReadsSD);
		}
	}//end of test
	
	public double TotalCounts(double[] array){	
		int totalCounts = 0;
		for (int i = 0; i <array.length;i++){totalCounts+=array[i];}
		return totalCounts;
	}
	
	public void printArray(double[] array){
		for (int i = 0; i < array.length; i++){System.out.print(array[i]+"\t");}
		System.out.println();
	}
	
	public double[] randomelyAssignReads(double[] composite){	
			
		double randomReads[] = new double[composite.length];
		for (int i = 0 ; i < randomReads.length; i++){ 
			randomReads[i] = 0;
		}	
		double totalCounts = TotalCounts(composite);	
		for (int i = 0 ; i < totalCounts; i++){
			randomReads[ThreadLocalRandom.current().nextInt(0, composite.length)] += 1;
		}				
		return randomReads;
	}
	
	public double[] shuffleComposite(double[] composite){
		
		double[] shuffledComposite =  new double[composite.length];
		System.arraycopy(composite, 0, shuffledComposite, 0, composite.length);
		
		Random rand = ThreadLocalRandom.current();
		for (int i = shuffledComposite.length - 1; i >0; i--){
			int index = rand.nextInt(i+1);
			double counts = shuffledComposite[index];
			shuffledComposite[index] = shuffledComposite[i];
			shuffledComposite[i] = counts;				
		}	
		return shuffledComposite;	 
	}
	
	// calculate weighted standard deviation
	public double[] computeWeightedStandardDeviation(double[] composite){
		
		double sumWeightedVar = 0 ; double sumWeights = 0; double weightedVar = 0; double weightedSD = 0 ;
		double [] distributionScore = new double [2]; 		
	
		int N = composite.length;
		double M = TotalCounts(composite);
		
		int maxIndex = -1;
		double maxCount = 0;
		for (int i = 0 ; i < N; i ++){
			if (composite[i] > maxCount){
				maxCount = composite[i];
				maxIndex = i;
			}
		}			
		for (int i = 0; i < N ; i++){
			sumWeightedVar += composite[i]*(i-maxIndex)*(i-maxIndex);
			sumWeights += composite[i];
		}
		weightedVar = (M*sumWeightedVar)/((M-1)*sumWeights);
		weightedSD = Math.sqrt(weightedVar);				
		distributionScore[0] = maxIndex - N/2;
		distributionScore[1] = weightedSD;
		
		return distributionScore;
	}
	
	// calculate standard deviation
	public double computeStandardDeviation(double[] x){
		
		int N = x.length;
		double var = 0 ; double sd = 0;	double sum = 0;	double mu = 0;
		for (int i = 0; i < N; i++){
			sum += x[i];
		}
		mu = sum/N;	
		for (int i = 0; i < N ; i++){
			var += (x[i]-mu)*(x[i]-mu);
		}
		sd = Math.sqrt(var/N);
		
		return sd;
	}
	
	public void printWeightedStandardDeviationStatistics(){		
		System.out.println("#sampleName\tMaxPos\tsampleWeightedSD\tp_val\tSignificant?\tz_score\tnullMu\tnullSD");
		for (Sample sample : sampleStandardDeviation.keySet()){			
			double [] distributionScore = sampleStandardDeviation.get(sample);		
			double [] stats = sampleStatsAndNull.get(sample);
			System.out.println(sample.getName()+"\t"+distributionScore[0]+"\t"+distributionScore[1]+"\t"+stats[0]+"\t"+stats[1]+"\t"+stats[2]+"\t"+stats[3]+"\t"+stats[4]);		
		}		
	}	
	
	public static void main(String[] args){
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		
		// parsing command line arguments
		ArgParser ap = new ArgParser(args);		
		int win = Args.parseInteger(args, "win", 1000);
		int fivePrimeShift = 0;
		fivePrimeShift = Args.parseInteger(args,"readshift", 6);
		List<StrandedPoint> spoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		FeatureCountsLoader fcLoader = new FeatureCountsLoader(gconf, econf, manager);
		fcLoader.setStrandedPoints(spoints);
		fcLoader.setWindowSize(win);
		fcLoader.setFivePrimeShift(fivePrimeShift);
		if (fivePrimeShift !=0){fcLoader.setEdge();}

		ReadDistributionProfiler profile = new ReadDistributionProfiler(fcLoader); 	
		profile.StandardDeviationFromExpAndNull();
		profile.StandardDeviationFromRandomReads();
		profile.printWeightedStandardDeviationStatistics();		
		if (Args.parseFlags(args).contains("printComposite")){profile.turnOnPrintComposite();}
		
		manager.close();
	}
}