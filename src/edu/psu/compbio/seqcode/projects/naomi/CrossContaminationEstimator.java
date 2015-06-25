package edu.psu.compbio.seqcode.projects.naomi;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;

public class CrossContaminationEstimator {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	
	public CrossContaminationEstimator(GenomeConfig gcon, ExptConfig econ){	
		gconfig = gcon;
		econfig = econ;
	}
	
	public float [][] getXYpairs(){
		
		ExperimentManager manager = new ExperimentManager(econfig);
		Genome genome = gconfig.getGenome();
	
		// a static array that can store all genome positions
		float[] [] dataPoints = new float [(int) genome.getGenomeLength()] [3];		
		for (int i = 0; i< genome.getGenomeLength(); i ++){
			for (int j = 0; j <3; j ++)
				dataPoints[i][j] = 0;		
		}
	
		int dataPointsIndex = 0;

		Iterator<Region> chroms = new ChromosomeGenerator<Genome>().execute(genome);
		//iterating each chromosome (each chromosome is a region).
		while (chroms.hasNext()) {
			
			Region currChrom = chroms.next();
			
			//get StrandedBaseCount list for each chromosome
			Map<Sample, List<StrandedBaseCount>> sampleCountsMap = new HashMap<Sample, List<StrandedBaseCount>>();
			for (Sample sample : manager.getSamples())
				sampleCountsMap.put(sample,sample.getBases(currChrom)); 			
			
			int currchromSize= currChrom.getWidth()+1;
			
			float[][] bpCounts = new float [currchromSize][sampleCountsMap.size()];
			for (int i = 0; i<currchromSize;i++){
				for (int j = 0; j<sampleCountsMap.size(); j++)
					bpCounts[i][j] = 0;			
			}

			//StrandedBasedCount object contains positive and negative strand separately
			for (Sample sample : manager.getSamples()){				
				List<StrandedBaseCount> currentCounts = sampleCountsMap.get(sample);
				for (StrandedBaseCount hits: currentCounts)
					bpCounts[hits.getCoordinate()][sample.getIndex()]+=hits.getCount();	
//					if (hits.getCount()>13) System.out.println("counts bigger than 13: "+hits.getCount());
				
				currentCounts = null;
			}

			int maxIndex = 0;
			float maxcounts = 0;
			float restSum = 0;
			
			//iterating bpCounts;  for each position in the chromosome, find the max among samples and copy it to dataPoints
			for (int i = 0; i<currchromSize;i++){
				// iterating one position among different samples
				for (int samp = 0; samp<sampleCountsMap.size();samp++){
					if (bpCounts[i][samp]>maxcounts){
						maxcounts = bpCounts[i][samp];
						maxIndex = samp;
					}					
				}
				
				for (int samp = 0; samp<sampleCountsMap.size();samp++){
					if (bpCounts[i][samp] != maxcounts)
						restSum = restSum + bpCounts[i][samp];					
				}
				
				if(maxcounts!=0){
					dataPoints[dataPointsIndex][0] = maxcounts;
					dataPoints[dataPointsIndex][1] = restSum;
					dataPoints[dataPointsIndex][2] = maxIndex;	
					
					dataPointsIndex++;
				}				
				maxIndex = 0;
				maxcounts = 0;
				restSum = 0;
			}
			bpCounts = null;			
			// remove all mapping 
			sampleCountsMap.clear();		
		}//end of chromosome iteration
		
		//iterate dataPoints to figure out non-zero dataPoints[i][0] (=maxTag number)
		int dataPointsSize = 0;		
		for (int i = 0; i<(int) genome.getGenomeLength(); i++){
			if (dataPoints[i][0]!=0)
				dataPointsSize++;
		}
		
		//copy non-zero dataPoints to xyPoints
		float [][] xyPairs = new float[dataPointsSize][3];		
		for (int i = 0; i<(int) genome.getGenomeLength();i++){
			if (dataPoints[i][0]!=0){
				for (int s= 0; s<3;s++)
					xyPairs[i][s]=dataPoints[i][s];
			}
		}
	
		return (xyPairs);
	}

	public void printXYpairs(){
			
		float [][] xyPairs = getXYpairs();
		
		//printing xyPairs
		System.out.println("#max_tag_number\tsum_of_other_sample's_tags\tsampleID");
		for (int i = 0; i< xyPairs.length; i++)
				System.out.println(xyPairs[i][0]+"\t"+xyPairs[i][1]+"\t"+xyPairs[i][2]);			
	}			
	
//	public void K_LineMeans(int k){
		
//		float [][] xyPairs = getXYpairs();
		
		
		
//	}
	
	public static void main(String[] args){
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		
		CrossContaminationEstimator estimator = new CrossContaminationEstimator (gconf, econf);
		estimator.getXYpairs();	
		estimator.printXYpairs();
		// trying k = 2
		int k = 2;		
//		estimator.K_LineMeans(k);
	}
}