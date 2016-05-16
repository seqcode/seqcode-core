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
	
	Map<Sample,double[]> sampleComposite = new HashMap<Sample,double[]>();
	Map<Sample,double[]> sampleStandardDeviation = new HashMap<Sample,double[]>();
	Map<Sample,Double> nullStandardDeviation = new HashMap<Sample,Double>();
	
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
	
	public void getCountsfromStrandedRegions(){
		
		// Because shift will shorten the array, edge will ensure that windows are covered with reads
		int edge = 40;
		
		// StrandedBaseCount list for each stranded regions for each sample
		Map<Sample, Map<StrandedRegion,List<StrandedBaseCount>>> sampleCountsMap = new HashMap<Sample, Map<StrandedRegion,List<StrandedBaseCount>>>();
		
		// sample counts array for each stranded region for each sample
		Map<Sample, Map<StrandedRegion,double[][]>> sampleCountsArray = new HashMap<Sample, Map<StrandedRegion,double[][]>>();
				
		List<StrandedRegion> regionList = new ArrayList<StrandedRegion>();
		for(Point p: strandedPoints){		
			int start = Math.max(1, p.getLocation() - window/2);
			int end = Math.min(p.getLocation() + window/2, p.getGenome().getChromLength(p.getChrom()));				
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
				
				System.out.println("size of regionCounts is "+regionCounts.length);
				for (int i = 0 ; i <50; i++){
					System.out.print(regionCounts[i][0]+"\t");
				}
				
				// get shifted composite for forward and reverse strands
				for (int j = 0 ; j <=window ; j++){
					nonstrandedComposite[j] += regionCounts[j-fivePrimeShift+edge/2][0]; 
					nonstrandedComposite[j] += regionCounts[j+fivePrimeShift+edge/2][1];
				}
				System.out.println("nonstrandedComposite");
				for (int i = 0 ; i <50; i++){
					System.out.print(nonstrandedComposite[i]);
				}
				
			}
			sampleComposite.put(sample, nonstrandedComposite);		
		}
		manager.close();
	}
		
	public void StandardDeviationFromExpAndNull(){
		
		int iterations = 1000;
			
		for (Sample sample : sampleComposite.keySet()){
			// calculate standard deviation from sample
			double[] composite = sampleComposite.get(sample);			
			sampleStandardDeviation.put(sample, computeWeightedStandardDeviation(composite));
		
			//shuffle the data points and calculate standard deviation from null	
			double [] arrNullSD = new double [iterations]; 
			for (int itr = 0 ; itr < iterations; itr++){
				
				double[] shuffledComposite =  new double[composite.length];
				System.arraycopy(composite, 0, shuffledComposite, 0, composite.length);	
				
				Random rand = ThreadLocalRandom.current();
				for (int i = shuffledComposite.length - 1; i >0; i--){
					int index = rand.nextInt(i+1);
					double counts = shuffledComposite[index];
					shuffledComposite[index] = shuffledComposite[i];
					shuffledComposite[i] = counts;				
				}		
				arrNullSD[itr] = computeWeightedStandardDeviation(shuffledComposite)[1];
			}
			double sum = 0;
			for (int j = 0 ; j < arrNullSD.length; j ++){
				sum += arrNullSD[j];
			}
			double mu = sum/arrNullSD.length;
			double sd = computeStandardDeviation(arrNullSD);		
			double x = sampleStandardDeviation.get(sample)[1];
			double z_score = (x - mu)/sd;
			double p_val = Erf.erfc(Math.abs(z_score)/Math.sqrt(2));
			
			nullStandardDeviation.put(sample, mu);
			
			if ((x < mu) && p_val <0.05){
				System.out.println("significant with p-value of "+p_val+": z score is "+z_score);
			}else{
				System.out.println("not significant with p-value of "+p_val+": z score is "+z_score);
			}	
		}
	}
	
	// calculate weighted standard deviation
	public double[] computeWeightedStandardDeviation(double[] composite){
		
		int N = composite.length;
		double sumWeightedVar = 0 ; double sumWeights = 0; double weightedVar = 0; double weightedSD = 0 ;
		double [] distributionScore = new double [2]; 		
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
		weightedVar = (N*sumWeightedVar)/((N-1)*sumWeights);
		weightedSD = Math.sqrt(weightedVar);				
		distributionScore[0] = N/2 - maxIndex;
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
	
	public void printSampleStandardDeviation(){		
		for (Sample sample : sampleStandardDeviation.keySet()){			
			double [] distributionScore = sampleStandardDeviation.get(sample);			
			System.out.println(sample.getName()+"\t"+distributionScore[0]+"\t"+distributionScore[1]);	
			double nullScore = nullStandardDeviation.get(sample);	
			System.out.println("null_score\t"+nullScore);	
		}		
	}	
	
	public void printSampleComposite(){		
		for (Sample sample : sampleComposite.keySet()){			
			double [] composite = sampleComposite.get(sample);
			System.out.println(sample.getName());
			for (int i = 0; i < composite.length; i++){
				System.out.print(composite[i]+"\t");
			}
			System.out.println();			
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
		profile.printSampleStandardDeviation();		
		if (Args.parseFlags(args).contains("printComposite")){profile.printSampleComposite();}
		
		manager.close();
	}
}