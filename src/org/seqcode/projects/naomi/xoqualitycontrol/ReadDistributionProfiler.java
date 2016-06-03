package org.seqcode.projects.naomi.xoqualitycontrol;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.special.Erf;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.projects.naomi.sequtils.FeatureCountsLoader;


/**
 * Holds methods to quantify non-uniform distributions around reference features
 * - Weighted standard deviation around the max counts as in Rohit thesis
 * - F statistics, which compares variance of control over samples
 * - Relative entropy measure for randomness
 * - Information content : difference between the entropy before and after the message I(X) = H_before - H_after
 * 
 * input : reference point
 * 
 * @author naomi yamada
 */

public class ReadDistributionProfiler {
	protected ExperimentManager manager;
	protected ExptConfig exptConfig;
	protected FeatureCountsLoader featureCountsLoader;
		
	protected List<StrandedPoint> strandedPoints;
	protected List<StrandedRegion> strandedRegions;
		
	protected int fivePrimeShift = 6;
	protected int window = 1000;
	protected int iterations = 1000;	
	protected boolean printSampleComposite = false;
	protected boolean zscore = false;

	protected Map<ControlledExperiment,double[]> sampleStandardDeviation = new HashMap<ControlledExperiment,double[]>();
	protected Map<Sample,double[]> sampleStatsAndNull = new HashMap<Sample,double[]>();
	protected Map<Sample,Double> Fstatistics = new HashMap<Sample,Double>();
	protected Map<Sample,Double> sampleRelativeEntropy = new HashMap<Sample,Double>();
	protected Map<Sample,double[]> informationContent = new HashMap<Sample,double[]>();
	
	public ReadDistributionProfiler(ExptConfig econf, ExperimentManager man, FeatureCountsLoader fcloader){	
		featureCountsLoader = fcloader;
		exptConfig = econf;
		manager = man;
	}
	
	// setters
	public void turnOnPrintComposite(){printSampleComposite = true;}
	public void turnOnZscore(){zscore = true;}
		
	public void computeReadDistributionStatistics(){
		
		double[] composite = new double[window+1];
		double[] contComposite = new double[window+1];
		Map<ControlledExperiment,double[]> controlStandardDeviation = new HashMap<ControlledExperiment,double[]>();
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){
				composite = featureCountsLoader.sampleComposite(rep); 
				contComposite = featureCountsLoader.controlComposite();
					
				// calculate standard deviation from sample
				if (printSampleComposite==true){ // print composites
					System.out.println(rep.getSignal().getName());
					printArray(composite);
					System.out.println(rep.getControl().getName());
					printArray(contComposite);
				}
			
				/// relative entropy calculation
				double [] fracSampleProfile = counts2Probability(composite);
				double [] fracContProfile = counts2Probability(contComposite);			
				sampleRelativeEntropy.put(rep.getSignal(), computeRelativeEntropy(fracSampleProfile,fracContProfile));	
			
				// information content calculation
				informationContent.put(rep.getSignal(), computeInformation(fracSampleProfile, fracContProfile));
			
				///normalizing counts
				double scaling = rep.getControlScaling();			
				double[] normalizedComposite = new double[composite.length];
				normalizedComposite[0] = composite[0] - scaling*(contComposite[0]+contComposite[1])*1/2;
				normalizedComposite[composite.length-1] = composite[composite.length-1] - scaling*(contComposite[0]+contComposite[1])*1/2;		
				for (int i = 1 ; i <composite.length-1 ; i++){
					normalizedComposite[i] = composite[i] - scaling*(contComposite[i-1]+contComposite[i]+contComposite[i+1])*1/3;
				}			
//				System.out.println("normalized composite ");
//				printArray(normalizedComposite);					
			
				sampleStandardDeviation.put(rep, computeWeightedStandardDeviation(composite));
				controlStandardDeviation.put(rep, computeWeightedStandardDeviation(contComposite));
						
				// F statistics
				double x = sampleStandardDeviation.get(rep)[1];
				double controlSD = controlStandardDeviation.get(rep)[1];
				double sampleVar = x*x;			
				double controlVar = controlSD*controlSD;		
				double Fstat = controlVar/sampleVar;
				Fstatistics.put(rep.getSignal(), Fstat);		
//				FDistribution fdist = new FDistribution(window-1, window-1);
//				double Fpval = 1- fdist.cumulativeProbability(Fstat);
//				System.out.println(rep.getSignal().getName()+": F statistics p-val is "+Fpval);		
				
				if (zscore == true){Zscore(rep,composite);}				
			}
		}
	}
	
	public void Zscore(ControlledExperiment rep, double[] composite){
		
		//shuffle the data points and calculate standard deviation from null	
		double [] arrNullSD = new double [iterations]; 
		for (int itr = 0 ; itr < iterations; itr++){				
			double[] shuffledComposite =  shuffleComposite(composite);		
			if (itr == 1){printArray(shuffledComposite);}
			arrNullSD[itr] = computeWeightedStandardDeviation(shuffledComposite)[1];
		}		
		//for test purpose
//		System.out.println("sd from null distribution");
//		printArray(arrNullSD);	
		
		double mu = TotalCounts(arrNullSD)/arrNullSD.length;
		double sd = computeStandardDeviation(arrNullSD);		
		double x = sampleStandardDeviation.get(rep)[1];		
		double minNullSD = 100000;
		for (int i = 0; i < arrNullSD.length; i++){
			if (arrNullSD[i] <minNullSD){
				minNullSD = arrNullSD[i];
			}
		}		
		/// Z score calculation
		double z_score = (x - mu)/sd;
		double p_val = Erf.erfc(Math.abs(z_score)/Math.sqrt(2));		
		double [] statistics = new double [5];
		statistics[0] = p_val;
		if ((x < mu) && p_val <0.05){statistics[1]=1;}else{statistics[1] = 0;}
		statistics[2] = z_score;
		statistics[3] = mu;
		statistics[4] = sd;					
		sampleStatsAndNull.put(rep.getSignal(), statistics);		
	}// End of Z score calculation	
	
	// This is for test purpose
	public void StandardDeviationFromRandomReads(){		
		for (ControlledExperiment rep : sampleStandardDeviation.keySet()){
			double[] composite = sampleStandardDeviation.get(rep);	
			double [] randomReadsSD = new double[iterations];
			for (int itr = 0 ; itr <iterations; itr++){
				double[] randomReads = randomelyAssignReads(composite);
				randomReadsSD[itr] = computeWeightedStandardDeviation(randomReads)[1];
			}
//			System.out.println("standard deviation from random reads");
//			printArray(randomReadsSD);
		}
	}//end of test
	
	/**************************
	**  Calculation Methods  **
	***************************/
	public double TotalCounts(double[] array){	
		int totalCounts = 0;
		for (int i = 0; i <array.length;i++){totalCounts+=array[i];}
		return totalCounts;
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
			sumWeightedVar += Math.abs(composite[i])*(i-maxIndex)*(i-maxIndex);
			sumWeights += Math.abs(composite[i]);
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
		double var = 0 ; double sd = 0;
		double mu = TotalCounts(x)/N;
		for (int i = 0; i < N ; i++){
			var += (x[i]-mu)*(x[i]-mu);
		}
		sd = Math.sqrt(var/N);		
		return sd;
	}
	
	public double[] counts2Probability(double [] array){		
		double totalCounts=TotalCounts(array);
		double[] probability = new double [array.length];
		for (int i = 0; i <array.length; i++){
			probability[i] = array[i]/totalCounts;
		}
		return probability;
	}
	
	public double computeRelativeEntropy(double[] p, double [] q){
		double relativeEntropy = 0;
		for (int i = 0; i< p.length; i++){
			if (p[i] != 0){relativeEntropy += p[i]*Math.log(p[i]/q[i]);}
		}	
		return relativeEntropy;		
	}
	
	// returns I(X) in index 0, Hbefore in index 1, Hafter in index 2
	public double[] computeInformation(double[] p, double [] q){
		double[] information = new double[3]; 
		double negHbefore = 0 ; // control distribution
		double negHafter = 0 ; // signal distribution
		for (int i = 0 ; i < p.length; i++){
			if (q[i] != 0){negHbefore += q[i]*Math.log(q[i]);}
			if (p[i] != 0){negHafter += p[i]*Math.log(p[i]);}
		}
		information[0] = - negHbefore + negHafter;
		information[1] = - negHbefore;
		information[2] = - negHafter;
		return information;
	}	
	
	/***********************
	**  Printing methods  **
	************************/	
	public void printArray(double[] array){
		for (int i = 0; i < array.length; i++){System.out.print(array[i]+"\t");}
		System.out.println();
	}
	
	public void printWeightedStandardDeviationStatistics(){		
		System.out.println("#sampleName\tMaxPos\tsampleWeightedSD\tp_val\tSignificant?\tz_score\tnullMu\tnullSD");
		for (ControlledExperiment rep : sampleStandardDeviation.keySet()){			
			double [] distributionScore = sampleStandardDeviation.get(rep);		
			double [] stats = sampleStatsAndNull.get(rep.getSignal());
			System.out.println(rep.getSignal().getName()+"\t"+distributionScore[0]+"\t"+distributionScore[1]+"\t"+stats[0]+"\t"+stats[1]+"\t"+stats[2]+"\t"+stats[3]+"\t"+stats[4]);		
		}		
	}	
	
	public void printStatistics(){	
		System.out.println("#sampleName\tFstatistics\trelativeEntropy\tinformationContent\tHbefore\tHafter");
		for (Sample rep : sampleRelativeEntropy.keySet())			
			System.out.println(rep.getName()+"\t"+Fstatistics.get(rep)+"\t"+sampleRelativeEntropy.get(rep)+"\t"+
		informationContent.get(rep)[0]+"\t"+informationContent.get(rep)[1]+"\t"+informationContent.get(rep)[2]); 
	}
	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);		
        if((!ap.hasKey("peaks") && !ap.hasKey("regions")) ) { 
            System.err.println("Usage:\n " +
                               "ReadDistributionProfiler " +
                               "--species <organism;genome> OR\n" +
                               "--geninfo <genome info> AND --seq <path to seqs>\n" +
                               "--peaks <file containing coordinates of peaks> \n" +
                               "--win <window of sequence to take around peaks> \n" +
                               "\nOPTIONS:\n" +
                               "--readshift \n" +
                               "--printComposite \n" +
                               "--zScore \n" +
                               "");
            System.exit(0);
        }
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		
		// parsing command line arguments	
		int win = Args.parseInteger(args, "win", 1000);
		int fivePrimeShift = 0;
		fivePrimeShift = Args.parseInteger(args,"readshift", 6);
		List<StrandedPoint> spoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		FeatureCountsLoader fcLoader = new FeatureCountsLoader(gconf, spoints, win);
		fcLoader.setFivePrimeShift(fivePrimeShift);

		ReadDistributionProfiler profile = new ReadDistributionProfiler(econf,manager, fcLoader); 	
		if (Args.parseFlags(args).contains("printComposite")){profile.turnOnPrintComposite();}
		profile.computeReadDistributionStatistics();
		if (Args.parseFlags(args).contains("zScore")){profile.turnOnZscore();}
		profile.printStatistics();

		manager.close();
	}
}