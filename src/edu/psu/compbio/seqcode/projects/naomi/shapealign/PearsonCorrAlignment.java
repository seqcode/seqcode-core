package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

/**
 * PearsonCorrAlignment : use Pearson correlation to get the best alignment
 * 
 * @author naomi yamada
 *
 */

public class PearsonCorrAlignment {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
		
	protected List<Point> points;
	protected List<Region> regions;
	protected int window;
	
	protected double error = 0;
	protected double totalNum = 0;
	
	static final double MINIMUM_VALUE = -10000;
	
	protected Map<Sample, Map<Region,double[][]>> countsArray = new HashMap<Sample,Map<Region,double[][]>>();
	
	public PearsonCorrAlignment(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
	}
	
	// setters
	public void setPoints(List<Point> p){points = p;}
	public void setRegions(List<Region> reg){regions = reg;} 
	public void setWidth(int w){window = w;}
	public void setCountsArray(Map<Sample, Map<Region,double[][]>> sampleCounts){countsArray = sampleCounts;}
	
	// prints error rate
	public void printErrorRate(){System.out.println("error is "+error+ " totalNum is "+ totalNum+ " error rate is "+ (error/totalNum));}
	
	public void loadData(){
		
		List<Region> region = new ArrayList<Region>();
		for(Point p: points){
			region.add(p.expand(window));
		}
		setRegions(region);

		//get StrandedBaseCount list for each regions per sample
		Map<Sample, Map<Region,List<StrandedBaseCount>>> sampleCountsMap = new HashMap<Sample, Map<Region,List<StrandedBaseCount>>>();
		Map<Sample, Map<Region,double[][]>> sampleCountsArray = new HashMap<Sample, Map<Region,double[][]>>();
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){				
				Map<Region,List<StrandedBaseCount>> regionCounts =  new HashMap<Region,List<StrandedBaseCount>>();				
				for (Region reg : regions){
					regionCounts.put(reg, rep.getSignal().getBases(reg));
				}
				sampleCountsMap.put(rep.getSignal(),regionCounts);
			}					
		}
		
		//StrandedBasedCount object contains positive and negative strand separately
		for (Sample sample : sampleCountsMap.keySet()){
			
			Map<Region,double[][]> regionCounts = new HashMap<Region,double[][]>();
			
			for (Region reg : sampleCountsMap.get(sample).keySet()){			
				double[][] sampleCounts = new double[window*2+1][2];
				for (int i = 0;i <= window*2;i++){
					for (int s = 0; s<2; s++)
						sampleCounts[i][s] = 0;
				}				
				for (StrandedBaseCount hits: sampleCountsMap.get(sample).get(reg)){
					if (hits.getStrand()=='+'){
						sampleCounts[hits.getCoordinate()-reg.getStart()][0] = hits.getCount();
					}else{
						sampleCounts[hits.getCoordinate()-reg.getStart()][1] = hits.getCount();
					}
				}
				regionCounts.put(reg, sampleCounts);
			}
			
			sampleCountsArray.put(sample, regionCounts);
		}
		setCountsArray(sampleCountsArray);
	}
	
	public void excutePearsonCorrAlign(){
		
		for (Sample sample : countsArray.keySet()){
			for (int i = 0; i <regions.size();i++){		
			
//			for (int i = 0; i <1 ; i++){	
				for (int j = i+1; j <regions.size();j++){
					System.out.println("region is "+regions.get(j).getLocationString());
					pearsonCorrelation(sample, regions.get(i), regions.get(j));		
				}
			}
		}				
	}
	
	public void pearsonCorrelation(Sample sample, Region regA, Region regB){
		
		//get midpoints
//		double regAmid = regA.getMidpoint().getLocation();
//		double regBmid = regB.getMidpoint().getLocation();
		
		//get counts
		double [][] regACounts = countsArray.get(sample).get(regA);
		double [][] regBCounts = countsArray.get(sample).get(regB);
		
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
			double r = pc.computePearsonCorr();
			boolean fromReverse = pc.isReverse();
			
			if (r > maxCorr ){
				maxCorr = r;
				reverse = fromReverse;
				offset = slidingW - window/2;
			}					
		}
		System.out.println("offsets is "+offset);
		System.out.println("correlation is "+maxCorr);

		double[][] alignedRegB = new double[window+1][2];
		double max_b = MINIMUM_VALUE;
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				if (regBCounts[window/2+offset+i][s] > max_b){max_b = regBCounts[window/2+offset+i][s];}
			}
		}
		
		System.out.println("max_b is "+max_b);

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
		
		// incrementing error 
		totalNum += 1;
		if (offset !=0){
			error += 1;
		}
	}
	
	public static void main(String[] args){
		
		ArgParser ap = new ArgParser(args);		
		int win = Args.parseInteger(args, "win", 200);
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);
		
		ExperimentManager manager = new ExperimentManager(econf);

		PearsonCorrAlignment profile = new PearsonCorrAlignment(gconf, econf, manager); 		
		
		List<Point> points = RegionFileUtilities.loadPeaksFromPeakFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		profile.setPoints(points);
		profile.setWidth(win);
		profile.loadData();
		profile.excutePearsonCorrAlign();
		profile.printErrorRate();	
	}
}
