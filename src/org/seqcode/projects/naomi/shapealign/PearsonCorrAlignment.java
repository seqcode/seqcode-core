package org.seqcode.projects.naomi.shapealign;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.projects.naomi.sequtils.FeatureCountsLoader;


/**
 * PearsonCorrAlignment : use Pearson correlation to get the best alignment
 * 
 * @author naomi yamada
 *
 */

public class PearsonCorrAlignment {
	
	protected FeatureCountsLoader fcloader;
	protected ExperimentManager manager;
		
	protected List<StrandedRegion> strandeRegions;
	protected int window;
	
	protected boolean isWeighted = false;
	protected double error = 0;
	protected double totalNum = 0;	
	static final double MINIMUM_VALUE = -10000;
	
	protected Map<ControlledExperiment, Map<StrandedRegion, double[][]>> strandedRegionSampleCounts = new HashMap<ControlledExperiment, Map<StrandedRegion,double[][]>>(); 	
	protected double [][] pairwiseCorrelation;
	protected double [] offsetArray;
	
	public PearsonCorrAlignment(FeatureCountsLoader featureCountsLoader,ExperimentManager man){	
		fcloader = featureCountsLoader;
		manager = man;
	}

	// setters
	public void setWidth(int w){window = w;}
	public void setWeighted(){ isWeighted = true;}
	
	// prints error rate
	public void printErrorRate(){System.out.println("error is "+error+ " totalNum is "+ totalNum+ " error rate is "+ (error/totalNum));}
	
	public void excutePearsonCorrAlign(){
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){			
				strandedRegionSampleCounts.put(rep, fcloader.strandedRegionSampleCounts(rep));
			}
		}
		strandeRegions = fcloader.getStrandedRegions();

		// initialize pairwiseDistance
		pairwiseCorrelation = new double [strandeRegions.size()][strandeRegions.size()];
		for (int i = 0 ; i< strandeRegions.size(); i++){
			for (int j = 0 ; j < strandeRegions.size(); j++)
				pairwiseCorrelation[i][j] = 1;
		}			
		//initiazlie offsetArray
		offsetArray = new double [window+1];
		for (int i = 0 ; i <= window ; i++)
			offsetArray[i] = 0;
		
		// loop through all possible pairwise region to do Pearson correlation
		for (ControlledExperiment cexpt : strandedRegionSampleCounts.keySet()){
			for (int i = 0; i <strandeRegions.size();i++){		
			
//			for (int i = 0; i <1 ; i++){	
				for (int j = i+1; j <strandeRegions.size();j++){
//					System.out.println("region is "+regions.get(j).getLocationString());
					double corr = pearsonCorrelation(cexpt, strandeRegions.get(i), strandeRegions.get(j));	
					pairwiseCorrelation[i][j] = corr;
					pairwiseCorrelation[j][i] = corr;
				}
			}
		}				
	}	
	
	public double pearsonCorrelation(ControlledExperiment controlledExpt, Region regA, Region regB){
		
		//get midpoints
//		double regAmid = regA.getMidpoint().getLocation();
//		double regBmid = regB.getMidpoint().getLocation();
		
		//get counts
		double [][] regACounts = strandedRegionSampleCounts.get(controlledExpt).get(regA);
		double [][] regBCounts = strandedRegionSampleCounts.get(controlledExpt).get(regB);
		
		//normalize the regAarray to set the max value 1
		// size of normRegACounts is [window_size][2]
		double maxA = MINIMUM_VALUE;
		for (int i = 0 ; i <=window ; i++){
			for (int s = 0 ; s <2 ; s++){
				if (regACounts[i+window/2][s] > maxA){maxA = regACounts[i+window/2][s];}
			}
		}
		double [][] normRegACounts = new double[window+1][2];
		for (int i = 0; i <=window; i++){
			for (int s = 0 ; s < 2 ; s++)
				normRegACounts[i][s] = regACounts[i+window/2][s]/maxA;
		}
		
		// compute Pearson correlation with sliding window
		double maxCorr = -1;
		boolean reverse = false;
		int offset = 0; 
		
		for (int slidingW = 0; slidingW <=window; slidingW++){			
			// copy and normalize B array
			double maxB = MINIMUM_VALUE;
			for (int i = 0; i <=window; i++){
				for (int s = 0 ; s < 2 ; s++){
					if (regBCounts[slidingW+i][s] > maxB){maxB = regBCounts[slidingW+i][s];}
				}			
			}
			double [][] normRegBCounts = new double [window+1][2];
			for (int i = 0; i <=window; i++){
				for (int s = 0 ; s < 2 ; s++){
					normRegBCounts[i][s] = regBCounts[slidingW+i][s]/maxB;
				}
			}
			
			//calculate Pearson Correlation
			PearsonCorrelation pc = new PearsonCorrelation(normRegACounts,normRegBCounts);
			
			double r = -1;
			boolean fromReverse = false;
			if (isWeighted !=true){
				r = pc.computePearsonCorr();
				fromReverse = pc.isReverse();
			}else{
				r = pc.computeWeightedPearsonCorr();
				fromReverse = pc.isReverse();				
			}	
			if (r > maxCorr ){
				maxCorr = r;
				reverse = fromReverse;
				offset = slidingW - window/2;
			}					
		}
//		System.out.println("offsets is "+offset);
//		System.out.println("correlation is "+maxCorr);

		double[][] alignedRegB = new double[window+1][2];
		double max_b = MINIMUM_VALUE;
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				if (regBCounts[window/2+offset+i][s] > max_b){max_b = regBCounts[window/2+offset+i][s];}
			}
		}		
//		System.out.println("max_b is "+max_b);

		for (int i = 0; i <= window; i++){
			for (int s = 0; s <2 ; s++){
				alignedRegB[i][s] = regBCounts[window/2+offset+i][s]/max_b;
			}
		}
		
		double temp_max_b = MINIMUM_VALUE;
		double [][] norm_b = new double[window+1][2];
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				if (regBCounts[window/2+i][s] > temp_max_b){temp_max_b = regBCounts[window/2+i][s];}
			}
		}
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				norm_b[i][s] = regBCounts[window/2+i][s]/temp_max_b;
			}
		}		
/**	
		if (Math.abs(offset) > 1){ // printing in case where alignment fails
		System.out.println("before alignment reg A");
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(normRegACounts[i][0]+",");
		System.out.println();
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(normRegACounts[i][1]+",");	
		System.out.println();
		
		System.out.println("before alignment reg B");
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(norm_b[i][0]+",");
		System.out.println();
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(norm_b[i][1]+",");	
		System.out.println();
		
		System.out.println("aligned reg B");
		if (reverse == true){
			for (int i = 0; i <alignedRegB.length; i++)
				System.out.print(alignedRegB[alignedRegB.length-1-i][0]+",");
			System.out.println();
			for (int i = 0; i <alignedRegB.length; i++)
				System.out.print(alignedRegB[alignedRegB.length-1-i][1]+",");						
		}else{
			for (int i = 0; i < alignedRegB.length; i++)
				System.out.print(alignedRegB[i][0]+",");
			System.out.println();
			for (int i = 0; i < alignedRegB.length; i++)
				System.out.print(alignedRegB[i][1]+",");			
		}		
		System.out.println();		
		}
**/		
		// increment offset array
		offsetArray[offset+window/2]++;				
		// incrementing error allowing offset of +-1
		totalNum += 1;
		if (Math.abs(offset) > 1){
			error += 1;
//			System.out.println(regA.getLocationString());
//			System.out.println(regB.getLocationString());
			}
		
		return maxCorr;	
	}
	
	public void printPairwiseCorrelation(){		
		for (int i = 0 ; i < strandeRegions.size(); i++){
			System.out.println(strandeRegions.get(i).getLocationString());
		}		
		for (int i = 0 ; i < strandeRegions.size();i++){
			for (int j = 0 ; j <strandeRegions.size(); j++){
				System.out.print(pairwiseCorrelation[i][j]+"\t");
			}
			System.out.println();	
		}
	}
	
	public void printOffsetArray(){
		System.out.println("offset array");
		for (int i = 0; i <= window ; i++)
			System.out.print(offsetArray[i]+"\t");
		System.out.println();
	}
	
	public static void main(String[] args){
				
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		
		// parsing command line arguments
		ArgParser ap = new ArgParser(args);		
		int win = Args.parseInteger(args, "win", 200);
		List<StrandedPoint> spoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		// window size must be twice bigger so it can slide window size times
		FeatureCountsLoader fcLoader = new FeatureCountsLoader(gconf, spoints, win*2);

		PearsonCorrAlignment profile = new PearsonCorrAlignment(fcLoader, manager); 	
		if (Args.parseFlags(args).contains("weighted")){profile.setWeighted();}
		profile.setWidth(win);
		profile.excutePearsonCorrAlign();
		profile.printErrorRate();	
		if (Args.parseFlags(args).contains("printCorrelation")){profile.printPairwiseCorrelation();}
		if (Args.parseFlags(args).contains("printOffset")){profile.printOffsetArray();}
		
		manager.close();
	}
}
