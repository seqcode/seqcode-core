package edu.psu.compbio.seqcode.projects.naomi.xoqualitycontrol;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.special.Erf;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

/**
 * ReadDistributionProfiler : calculates weighted standard deviation around the reference point as in Rohit thesis
 * 
 * input : reference point
 * 
 * @author naomi yamada
 */

public class ReadDistributionProfiler {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
		
	protected List<StrandedPoint> strandedPoints;
	protected List<StrandedRegion> strandedRegions;
	
	protected int fivePrimeShift = 6;
	protected int window = 1000;
	int iterations = 1000;	
	protected boolean printSampleComposite = false;
	
	Map<Sample,double[]> sampleComposite = new HashMap<Sample,double[]>();
	Map<Sample,double[]> sampleStandardDeviation = new HashMap<Sample,double[]>();
	Map<Sample,double[]> sampleStatsAndNull = new HashMap<Sample,double[]>();
	
	public ReadDistributionProfiler(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
	}
	
	// setters
	public void setStrandedPoints(List<StrandedPoint> p){strandedPoints = p;}
	public void setStrandedRegions(List<StrandedRegion> reg){strandedRegions = reg;} 
	public void setWindowSize(int w){window = w;}
	public void setFivePrimeShift(int s){fivePrimeShift = s;}
	public void turnOnPrintComposite(){printSampleComposite = true;}
	
	public void getCountsfromStrandedRegions(){
		
		// Because shift will shorten the array, edge will ensure that windows are covered with reads
		int edge = 40;
		
		// StrandedBaseCount list for each stranded regions for each sample
		Map<Sample, Map<StrandedRegion,List<StrandedBaseCount>>> sampleCountsMap = new HashMap<Sample, Map<StrandedRegion,List<StrandedBaseCount>>>();
		
		// sample counts array for each stranded region for each sample
		Map<Sample, Map<StrandedRegion,double[][]>> sampleCountsArray = new HashMap<Sample, Map<StrandedRegion,double[][]>>();
				
		List<StrandedRegion> regionList = new ArrayList<StrandedRegion>();
		for(Point p: strandedPoints){		
			int start = Math.max(1, p.getLocation() - (window+edge)/2 );
			int end = Math.min(p.getLocation() + (window+edge)/2, p.getGenome().getChromLength(p.getChrom()));				
			StrandedRegion strandedReg = new StrandedRegion(p.getGenome(), p.getChrom(), start, end, p.getStrand());					
			regionList.add(strandedReg);
		}
		setStrandedRegions(regionList);
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){				
				Map<StrandedRegion,List<StrandedBaseCount>> regionCounts =  new HashMap<StrandedRegion,List<StrandedBaseCount>>();				
				for (StrandedRegion reg : strandedRegions){
					regionCounts.put(reg, rep.getSignal().getBases(reg));
				}
				sampleCountsMap.put(rep.getSignal(),regionCounts);
			}					
		}
		
		//StrandedBasedCount object contains positive and negative strand separately
		// Reverse the array depending of strand of features	
		for (Sample sample : sampleCountsMap.keySet()){
			
			Map<StrandedRegion,double[][]> regionCounts = new HashMap<StrandedRegion,double[][]>();
			
			for (StrandedRegion reg : sampleCountsMap.get(sample).keySet()){			
				double[][] sampleCounts = new double[window+edge+1][2];
				for (int i = 0;i <= window+edge;i++){
					for (int s = 0; s<2; s++)
						sampleCounts[i][s] = 0;
				}	
				if (reg.getStrand() == '+'){ // regions(features) are positive strand
					reg.getMidpoint().getLocation();
					
					for (StrandedBaseCount hits: sampleCountsMap.get(sample).get(reg)){	
						if (hits.getStrand()=='+'){
							sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][0] = hits.getCount();
						}else{
							sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][1] = hits.getCount();
						}					
					}
				}else{ // if regions (features) are reverse strand, I need to flip the strands and locations
					for (StrandedBaseCount hits: sampleCountsMap.get(sample).get(reg)){	
						if (hits.getStrand()=='+'){
							sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][1] = hits.getCount();
						}else{
							sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][0] = hits.getCount();	
						}			
					}		
				}
				regionCounts.put(reg, sampleCounts);
			}
			sampleCountsArray.put(sample, regionCounts);
		}
					
		// nonstrandedComposite is a shifted and strand merged composite to be used to measure standard deviation
		for (Sample sample : sampleCountsArray.keySet()){
			
			double [] nonstrandedComposite = new double[window+1];
			for (int i = 0; i <=window ; i ++)
				nonstrandedComposite[i] = 0;
			
			for (Region reg : sampleCountsArray.get(sample).keySet()){				
				double[][] regionCounts = sampleCountsArray.get(sample).get(reg);		

				// get shifted composite for forward and reverse strands
				for (int j = 0 ; j <=window ; j++){
					nonstrandedComposite[j] += regionCounts[j-fivePrimeShift+edge/2][0]; 
					nonstrandedComposite[j] += regionCounts[j+fivePrimeShift+edge/2][1];
				}
			}
			sampleComposite.put(sample, nonstrandedComposite);		
		}
		manager.close();
	}
		
	public void StandardDeviationFromExpAndNull(){

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
		int fivePrimeShift = Args.parseInteger(args,"readshift", 6);
		List<StrandedPoint> spoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);

		ReadDistributionProfiler profile = new ReadDistributionProfiler(gconf, econf, manager); 	
		profile.setStrandedPoints(spoints);
		profile.setWindowSize(win);
		profile.setFivePrimeShift(fivePrimeShift);
		profile.getCountsfromStrandedRegions();
		profile.StandardDeviationFromExpAndNull();
		profile.StandardDeviationFromRandomReads();
		profile.printWeightedStandardDeviationStatistics();		
		if (Args.parseFlags(args).contains("printComposite")){profile.turnOnPrintComposite();}
		
		manager.close();
	}
}