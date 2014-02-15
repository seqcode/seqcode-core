package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.deepseq.PairedCountData;
import edu.psu.compbio.seqcode.gse.utils.models.Model;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataFrame;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataRegression;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;

/**
 * ExperimentScaler: calculate a scaling transformation between two deep-seq experiments
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentScaler {
	protected Genome genome;
	protected Sample exptA, exptB;
	protected int windowSize=10000;
	protected double scalingRatio=-1;
	
	public ExperimentScaler(Sample a, Sample b){
		exptA = a;
		exptB = b;
		genome = a.getGenome();
	}
	
	/**
	 * Calculate a scaling ratio by fitting a line through the hit count pairs
	 * @return double 
	 */
	public double scalingRatioByRegression(int win){
		scalingRatio=1;
		if(exptB==null)
			return(1);
		windowSize = win;
		ArrayList<PairedCountData> scalingData = new ArrayList<PairedCountData>();
		for(String chrom:genome.getChromList()) {
            int chrlen = genome.getChromLength(chrom);
            for (int start = 1; start  < chrlen - windowSize; start += windowSize) {
                Region r = new Region(genome, chrom, start, start + windowSize);
                double countA = exptA.countHits(r);
                double countB = exptB.countHits(r);
                scalingData.add(new PairedCountData(countA , countB));                
            }
        }
		//Scaling ratio via Tim's regression                                                                                                                               
        DataFrame df = new DataFrame(edu.psu.compbio.seqcode.gse.deepseq.PairedCountData.class, scalingData.iterator());                                                                      
        DataRegression r = new DataRegression(df, "x~y - 1");                                                                                                               
        r.calculate();                                                                                                                                                      
        Map<String, Double> map = r.collectCoefficients();                                                                                                                  
        scalingRatio = map.get("y");                                                                                                                                        
        System.err.println(String.format("Scaling ratio estimated by regression = %.3f based on %d regions of size %d",
        		scalingRatio, scalingData.size(), windowSize));
        return(scalingRatio);
	}
	
	/**
	 * Find the median hit count ratio 
	 * @return
	 */
	public double scalingRatioByMedian(int win){
		scalingRatio=1;
		if(exptB==null)
			return(1);
		windowSize = win;
	    ArrayList<Float> ratios = new ArrayList<Float>();
	    //System.err.println("SCALING: "+exptA.getName()+" vs "+exptB.getName());
		for(String chrom:genome.getChromList()) {
            int chrlen = genome.getChromLength(chrom);
            for (int start = 0; start  < chrlen - windowSize; start += windowSize) {
                Region r = new Region(genome, chrom, start, start + windowSize);
                double countA = exptA.countHits(r);
                double countB = exptB.countHits(r);
                if(countB>0)
                	ratios.add((float)(countA / countB));
                else
                	ratios.add((float)(countA / 1));
                double tmpB = countB>0? countB:1;
                //System.err.println(countA+"\t"+tmpB);
            }
        }
        Collections.sort(ratios);
		scalingRatio = ratios.get(ratios.size() / 2);
        System.err.println(String.format("Scaling ratio estimated by median scaling = %.3f based on %d regions of size %d",
        		scalingRatio, ratios.size(), windowSize));
		
		return(scalingRatio);
	}
	
	/**
	 * Find the scaling ratio according to the SES method from Diaz, et al. Stat Appl Genet Mol Biol. 2012.
	 * Also sets a background proportion estimate for the signal channel.  
	 * @return
	 */
	public double scalingRatioBySES(int win){
		windowSize = win;
		scalingRatio=1;
		ArrayList<PairedCounts> counts = new ArrayList<PairedCounts>();
		
		//Either A alone if no experiment B, or A & B otherwise
		double totalA=0, totalB=0, numWin=0;
		for(String chrom:genome.getChromList()) {
            int chrlen = genome.getChromLength(chrom);
            for (int start = 0; start  < chrlen - windowSize; start += windowSize) {
                Region r = new Region(genome, chrom, start, start + windowSize);
                double countA = exptA.countHits(r);
                double countB = exptB==null ? 0 : exptB.countHits(r);                
                counts.add(new PairedCounts(countA, countB));
                
                totalA+=countA;
                totalB+=countB;
                numWin++;
            }
        }
		if(exptB==null){
			int winCount=0;
			totalB=totalA;
			for(String chrom:genome.getChromList()) {
	            int chrlen = genome.getChromLength(chrom);
	            for (int start = 0; start  < chrlen - windowSize; start += windowSize) {                
	                counts.get(winCount).y=totalA/numWin;
	                winCount++;
	            }
	        }
		}
		
        Collections.sort(counts);
        
        //SES procedure
        double cumulA=0, cumulB=0, maxDiffAB=0, maxDiffAprop=0, currDiff=0;
        int maxDiffIndex=0, i=0;
        for(PairedCounts pc : counts){
        	cumulA+=pc.x;
        	cumulB+=pc.y;
        	currDiff = (cumulB/totalB)-(cumulA/totalA);
        	if(currDiff>maxDiffAB){
        		maxDiffAB=currDiff;
        		maxDiffIndex=i;
        		maxDiffAprop=(cumulA/totalA);
        		scalingRatio = cumulB>0 ? cumulA/cumulB : 1;
        	}
        	i++;
        }
        
        System.err.println(String.format("Scaling ratio estimated by SES = %.3f based on %d regions of size %d",
        		scalingRatio, counts.size(), windowSize));
        
		return(scalingRatio);
	}
	
	/**
	 * Calculate the background proportion of an IP experiment by correcting the scaling ratio by the read count ratio.
	 * Be careful with this method, there are a couple of assumptions:
	 *  - The scaling ratio was calculated between IP and control experiments
	 *  - The method used to calculate the scaling ratio attempted to normalize to background and not all regions (e.g. SES method attempts background normalization) 
	 * @return
	 */
	public Double calculateBackgroundFromScalingRatio(){
		if(scalingRatio==-1)
			return(-1.0); //scaling not yet performed
		double ctrlCount = exptB==null ? exptA.getHitCount() : exptB.getHitCount();
		return(scalingRatio / (exptA.getHitCount()/ctrlCount));
	}
	
	/**
	 * Main for testing
	 * @param args
	 */
	public static void main(String[] args){
		
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("ExperimentScaler:");
			System.err.println(config.getArgsList());
		}else{
			Genome gen = config.getGenome();
			ExperimentManager exptMan = new ExperimentManager(config);
			ExperimentSet eset = exptMan.getExperimentSet();
			
			//Test
			System.err.println("Conditions:\t"+eset.getConditions().size());
			for(ExperimentCondition c : eset.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : eset.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			
			//Inter-sample normalization
			for(Sample sampA : eset.getSamples()){
				for(Sample sampB : eset.getSamples()){
					if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex()){
						
						System.out.println("\n"+sampA.getName()+" vs "+sampB.getName());
						double hitRatio = sampA.getHitCount()/sampB.getHitCount();
						System.out.println("Hit ratio ="+hitRatio);
						
						ExperimentScaler scaler = new ExperimentScaler(sampA, sampB);
						double medRatio = scaler.scalingRatioByMedian(10000);
						System.out.println("Scaling ratio by median = "+medRatio);
						
						double regRatio = scaler.scalingRatioByRegression(10000);
						System.out.println("Scaling ratio by regression = "+regRatio);
						
						double sesRatio = scaler.scalingRatioBySES(200);
						System.out.println("Scaling ratio by SES = "+sesRatio);
						
						double backProp = scaler.calculateBackgroundFromScalingRatio();
						System.out.println("Estimated background proportion = "+backProp);
					}
				}
			}			
		}
	}
	
	
	/**
	 * Simple class for storing paired counts that are sortable in first dimension
	 * @author mahony
	 *
	 */
	public class PairedCounts  implements Comparable<PairedCounts>{
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
