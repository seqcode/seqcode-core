package edu.psu.compbio.seqcode.deepseq.experiments;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.utils.models.Model;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataFrame;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataRegression;

/**
 * ExperimentScaler: calculate a scaling transformation between all Sample pairs in an ExperimentCondition
 * This is performed on the condition level so that a scaling can also be defined between pooled hits from all signals & controls
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentScaler {
	
	public ExperimentScaler(){}
	
	/**
	 * Calculate a scaling ratio by fitting a line through the hit count pairs.
	 * Using a 10Kbp window, this is the same as PeakSeq with Pf=0
	 * @return double 
	 */
	public double scalingRatioByRegression(List<Float> setA, List<Float> setB){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
			
		List<PairedCounts> scalingData = new ArrayList<PairedCounts>();
		for(int x=0; x<setA.size(); x++)
			scalingData.add(new PairedCounts(setA.get(x), setB.get(x)));                

		//Scaling ratio via Tim's regression                                                                                                                               
        DataFrame df = new DataFrame(PairedCounts.class, scalingData.iterator());                                                                      
        DataRegression r = new DataRegression(df, "x~y - 1");                                                                                                               
        r.calculate();                                                                                                                                                      
        Map<String, Double> map = r.collectCoefficients();                                                                                                                  
        scalingRatio = map.get("y");                                                                                                                                        
        return(scalingRatio);
	}
	
	/**
	 * Find the median hit count ratio 
	 * @return
	 */
	public double scalingRatioByMedian(List<Float> setA, List<Float> setB){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
			
		ArrayList<Float> ratios = new ArrayList<Float>();
	    for(int x=0; x<setA.size(); x++){
			if(setA.get(x)>0)
				ratios.add((float)(setA.get(x) / setB.get(x)));
			else
				ratios.add((float)(setA.get(x) / 1));
        }
        Collections.sort(ratios);
		scalingRatio = ratios.get(ratios.size() / 2);
        return(scalingRatio);
	}
	
	/**
	 * Find the scaling ratio according to the SES method from Diaz, et al. Stat Appl Genet Mol Biol. 2012.
	 * Also sets a background proportion estimate for the signal channel.  
	 * @return
	 */
	public double scalingRatioBySES(List<Float> setA, List<Float> setB){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
		
		float totalA=0, totalB=0;
		List<PairedCounts> counts = new ArrayList<PairedCounts>();
		for(int x=0; x<setA.size(); x++){
			totalA += setA.get(x);
			totalB += setB.get(x);
			counts.add(new PairedCounts(setA.get(x), setB.get(x)));                
		}
		
		Collections.sort(counts);
        
        //SES procedure
        double readRatio = totalA/totalB;
        double cumulA=0, cumulB=0, maxDiffAB=0, maxDiffAprop=0, currDiff=0;
        int maxDiffIndex=0, i=0;
        for(PairedCounts pc : counts){
        	cumulA+=pc.x;
        	cumulB+=pc.y;
        	currDiff = (cumulB/totalB)-(cumulA/totalA);
        	if(currDiff>maxDiffAB && cumulA>0 && cumulB>0){
        		maxDiffAB=currDiff;
        		maxDiffIndex=i;
        		maxDiffAprop=(cumulA/totalA);
        		scalingRatio = cumulA/cumulB;
        	}
        	i++;
        }
		return(scalingRatio);
	}
	
	/**
	 * Calculate the background proportion of an IP experiment by correcting the scaling ratio by the read count ratio.
	 * Be careful with this method, there are a couple of assumptions:
	 *  - The scaling ratio was calculated between IP and control experiments
	 *  - The method used to calculate the scaling ratio attempted to normalize to background and not all regions (e.g. SES method attempts background normalization) 
	 * @return
	 */
	public Double calculateBackgroundFromScalingRatio(ControlledExperiment expt){
		if(expt.getControlScaling()==-1)
			return(-1.0); //scaling not yet performed
		double ctrlCount = expt.getControl()==null ? expt.getSignal().getHitCount() : expt.getControl().getHitCount();
		return(expt.getControlScaling() / (expt.getSignal().getHitCount()/ctrlCount));
	}
	
	/**
	 * Main for testing
	 * @param args
	 */
	public static void main(String[] args){
		
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		
		if(gconfig.helpWanted()){
			System.err.println("ExperimentScaler:");
			System.err.println(gconfig.getArgsList()+"\n"+econfig.getArgsList());
		}else{
			Genome gen = gconfig.getGenome();
			ExperimentManager exptMan = new ExperimentManager(econfig);
			
			//Test
			System.err.println("Conditions:\t"+exptMan.getConditions().size());
			for(ExperimentCondition c : exptMan.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : exptMan.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			
			ExperimentScaler scaler = new ExperimentScaler();
			
			//Generate the data structures for calculating scaling factors
			//10000bp window (median & regression)
			Genome genome = econfig.getGenome();
			Map<Sample, List<Float>> sampleWindowCounts = new HashMap<Sample, List<Float>>();
			for(Sample samp : exptMan.getSamples()){
				List<Float> currSampCounts = new ArrayList<Float>();
				for(String chrom:genome.getChromList()) {
		            int chrlen = genome.getChromLength(chrom);
		            for (int start = 1; start  < chrlen - 10000; start += 10000) {
		                Region r = new Region(genome, chrom, start, start + 10000);
		                currSampCounts.add(samp.countHits(r));
		            }
		        }
				sampleWindowCounts.put(samp, currSampCounts);
			}
			//Hit ratios
			for(Sample sampA : exptMan.getSamples())
				for(Sample sampB : exptMan.getSamples())
					if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex()){
						double hitRatio = sampA.getHitCount()/sampB.getHitCount();
						System.out.println("HitRatio\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+hitRatio);
					}
			
			//Median
			for(Sample sampA : exptMan.getSamples())
				for(Sample sampB : exptMan.getSamples())
					if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
						System.out.println("Median\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByMedian(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB)));
			//Regression
			for(Sample sampA : exptMan.getSamples())
				for(Sample sampB : exptMan.getSamples())
					if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
						System.out.println("Regression\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByRegression(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB)));

			
			//1000bp window (SES)
			sampleWindowCounts = new HashMap<Sample, List<Float>>();
			for(Sample samp : exptMan.getSamples()){
				List<Float> currSampCounts = new ArrayList<Float>();
				for(String chrom:genome.getChromList()) {
		            int chrlen = genome.getChromLength(chrom);
		            for (int start = 1; start  < chrlen - 1000; start += 1000) {
		                Region r = new Region(genome, chrom, start, start + 1000);
		                currSampCounts.add(samp.countHits(r));
		            }
		        }
				sampleWindowCounts.put(samp, currSampCounts);
			}
			//SES
			for(Sample sampA : exptMan.getSamples())
				for(Sample sampB : exptMan.getSamples())
					if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
						System.out.println("SES\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioBySES(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB)));

			exptMan.close();
		}
	}
	
	
	/**
	 * Simple class for storing paired counts that are sortable in first dimension
	 * @author mahony
	 *
	 */
	public class PairedCounts extends Model implements Comparable<PairedCounts>{
		public Double x,y;
		public PairedCounts(double a, double b){
			x=a;
			y=b;
		}
		//Sort on the X variables
		public int compareTo(PairedCounts pc) {
			if(x<pc.x){return -1;}
			if(x>pc.x){return 1;}
			return 0;
		}		
	}
}
