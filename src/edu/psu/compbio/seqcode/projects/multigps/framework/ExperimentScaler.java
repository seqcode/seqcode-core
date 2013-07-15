package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.deepseq.PairedCountData;
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
		double scalingRatio=1;
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
        System.err.println(String.format("Scaling ratio estimated as %.3f based on %d regions of size %d",
        		scalingRatio, scalingData.size(), windowSize));
        return(scalingRatio);
	}
	
	/**
	 * Find the median hit count ratio 
	 * @return
	 */
	public double scalingRatioByMedian(int win){
		double scalingRatio=1;
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
        System.err.println(String.format("Scaling ratio estimated as %.3f based on %d regions of size %d",
        		scalingRatio, ratios.size(), windowSize));
		
		return(scalingRatio);
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
					}
				}
			}			
		}
	}
}
