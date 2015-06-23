package edu.psu.compbio.seqcode.projects.naomi;

import java.util.ArrayList;
//import java.util.Collections;
//import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;
//import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class CrossContaminationEstimator {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
//	protected ExperimentManager manager;
	
//	public CrossContaminationEstimator(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){
	public CrossContaminationEstimator(GenomeConfig gcon, ExptConfig econ){	
		gconfig = gcon;
		econfig = econ;
//		manager = man;
	}
	
	public void printDataPoints(){
		
		ExperimentManager manager = new ExperimentManager(econfig);

		Genome genome = gconfig.getGenome();	
//		List<String> chromNames = genome.getChromList();	

		int sampleSize = manager.getSamples().size();
		
		float[] tempCounts = new float [sampleSize];
		
		// this is the static array that can store all genome positions
		float[] [] dataPoints = new float [(int) genome.getGenomeLength()] [3];
		
		//initialize dataPoints array
		for (int i = 0; i< genome.getGenomeLength(); i ++){
			for (int j = 0; j <3; j ++){
				dataPoints[i][j] = 0;
			}
		}
		
		int dataPointsIndex = 0;
		
		Iterator<Region> chroms = new ChromosomeGenerator<Genome>().execute(genome);
		//iterating each chromosome; each chromosome is a region.
		while (chroms.hasNext()) {
			
			Region currChrom = chroms.next();
			
			//get base pair counts for each chromosome
			List<StrandedBaseCount> bpcounts = new ArrayList<StrandedBaseCount>();
			for (Sample sample : manager.getSamples()){
				bpcounts = sample.getBases(currChrom); 
			}
			
			int maxIndex = 0;
			float maxcounts = 0;
			float restSum = 0;

			for (int i = 0; i < bpcounts.size(); i ++){
				
				for (int j = 0; sampleSize <0; j++){
					tempCounts[j] = bpcounts.get(i).getCount();
				}			
			
				for (int j = 0; sampleSize <0; j++){
					if (tempCounts[j] > maxcounts){
						maxcounts = tempCounts[j];
						maxIndex = j;
					}
				}
				for (int j = 0; sampleSize <0; j++){
					if (maxcounts != tempCounts[j]){
						restSum = restSum + tempCounts[j];
					}
				}
				if (maxcounts!= 0){
					dataPoints[dataPointsIndex][0] = maxcounts;
					dataPoints[dataPointsIndex][1] = restSum;
					dataPoints[dataPointsIndex][2] = maxIndex;
				}
				dataPointsIndex++;
				maxcounts = 0;
				restSum = 0;		
				
				for (int j = 0; sampleSize<0;j++){
					tempCounts[j]= 0;
				}
			}
		}
		
		//printing datapoints
		int i = 0;
		while (dataPoints [i][0] !=0){
			System.out.println(dataPoints[i][0]+"\t"+dataPoints[i][1]+"\t"+dataPoints[i][2]);			
			}
			i ++;
	}							
	
	
	public static void main(String[] args){
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		
//		ArgParser ap = new ArgParser(args);
		
		CrossContaminationEstimator estimator = new CrossContaminationEstimator (gconf, econf);
		estimator.printDataPoints();
		
	}
}
