package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
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

public class SmithWatermanAlignment {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	
	protected List<Point> points;
	protected List<Region> regions;
	protected int window;
	
	//constants for Smith-Waterman Algorithms
	final static float GAP_OPEN = 2;
	final static float GAP_EXT = 1;
	final static float MINIMUM_VALUE = Float.MIN_VALUE;
	
	protected Map<Sample, Map<Region,float[][]>> countsArray;
	
	public SmithWatermanAlignment(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
	}
	
	public void setPoints(List<Point> p){points = p;}
	public void setRegions(List<Region> reg){regions = reg;} 
	public void setWidth(int w){window = w;}
	
	public float computeScore(float aVal, float bVal){
		float score = (aVal + bVal)/2 - Math.abs(aVal - bVal);
		return score;		
	}
	
	public void smithWatermanAlgorithm(Sample sample, Region regA, Region regB){
		
		//get counts
		float [][] regACounts = countsArray.get(sample).get(regA);
		float [][] regBCounts = countsArray.get(sample).get(regB);
		
		//normalize the arrays to set the max value 1
		//should I be normalizing using max of either strand ?
		float maxA = MINIMUM_VALUE;
		float maxB = MINIMUM_VALUE;
		for (int i = 0; i <window ; i++){
			for (int s = 0 ; s < 2 ; s++){
				if (regACounts[i][s] > maxA){maxA = regACounts[i][s];}
				if (regBCounts[i][s] > maxB){maxB = regBCounts[i][s];}
			}
		}
		
		float [][] normRegACounts = new float [window][2];
		float [][] normRegBCounts = new float [window][2];
		for (int i = 0 ; i <window; i++){
			for (int s = 0 ; s < 2 ; s++){
				normRegACounts[i][s] = regACounts[i][s]/maxA;
				normRegBCounts[i][s] = regBCounts[i][s]/maxB;
			}
		}
				
		//initialization of M[][] & I[][] matrix
		float [][] M = new float [window+1][window+1];
		float [][] I = new float [window+1][window+1];
		char [][] H = new char [window+1][window+1];
		for (int i = 0 ; i < window; i++){
			for (int j = 0 ; j < window; j++){
				M[i][j] = 0;
				I[i][j] = 0;
				H[i][j] = '.';
			}
		}
		
		//fill in M[i][j] & I[i][j] matrix
		for (int i = 1 ; i <window; i++){
			for (int j = 1 ; j <window ; j++){
				float mScore = computeScore(normRegACounts[i-1][0], normRegBCounts[i-1][0])
						+ computeScore(normRegACounts[i-1][1], normRegBCounts[i-1][1]);
				
				M[i][j] = Math.max(M[i-1][j-1] + mScore, I[i-1][j-1] + mScore);
				
				float temp_I[] = new float[3];
				temp_I[0] = M[i][j-1] - GAP_OPEN;
				temp_I[1] = I[i][j-1] - GAP_EXT;
				temp_I[2] = M[i-1][j] - GAP_OPEN;
				temp_I[3] = I[i-1][j] - GAP_EXT;
				
				float max_I = MINIMUM_VALUE;
				for (int k = 0 ; i < 4 ; i++){
					if (temp_I[k] > max_I){ max_I = temp_I[i];}
				}
				
				I[i][j] = max_I;
			}
			
			// find the highest value
			
			
			// back track to reconstruct the path
				
		}
		
		
	}
	
	public void loadData(){
		
		List<Region> regions = new ArrayList<Region>();
		for(Point p: points){
			regions.add(p.expand(window/2));
		}
		
		//merge overlapping regions and set back to the original size				
		List<Region> resizedRegions = new ArrayList<Region>();		
		for (Region reg : Region.mergeRegions(regions)){ //mergeRegions is a static method	
			resizedRegions.add(reg.resize(window));
		}		
		setRegions(resizedRegions);
		
		//get StrandedBaseCount list for each regions per sample
		Map<Sample, Map<Region,List<StrandedBaseCount>>> sampleCountsMap = new HashMap<Sample, Map<Region,List<StrandedBaseCount>>>();
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){				
				Map<Region,List<StrandedBaseCount>> regionCounts =  new HashMap<Region,List<StrandedBaseCount>>();				
				for (Region reg : resizedRegions){
					regionCounts.put(reg, rep.getSignal().getBases(reg));
				}
				sampleCountsMap.put(rep.getSignal(),regionCounts);
			}					
		}
		
		//StrandedBasedCount object contains positive and negative strand separately
		for (Sample sample : sampleCountsMap.keySet()){
			
			Map<Region,float[][]> regionCounts = new HashMap<Region,float[][]>();
			
			for (Region reg : sampleCountsMap.get(sample).keySet()){			
				float[][] sampleCounts = new float[window][2];
				for (int i = 0;i<window;i++){
					for (int s = 0; s<2; s++)
						sampleCounts[i][s] = 0;
				}				
				for (StrandedBaseCount hits: sampleCountsMap.get(sample).get(reg)){
					if (hits.getStrand()=='+'){
						sampleCounts[hits.getCoordinate()][1] = hits.getCount();
					}else{
						sampleCounts[hits.getCoordinate()][2] = hits.getCount();
					}
				}
				regionCounts.put(reg, sampleCounts);
			}
			countsArray.put(sample, regionCounts);
		}
	}
	
	public void excuteShapeAlign(){
		
//		SmithWatermanAlignment pairAlign;
		for (Sample sample : countsArray.keySet()){
			for (int i = 0; i <regions.size();i++){			
				SmithWatermanAlignment pairAlign = null;
				for (int j = i+1; j <regions.size();j++)
					pairAlign.smithWatermanAlgorithm(sample, regions.get(i), regions.get(j));
			}
		}				
	}

	
	public static void main(String[] args){
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		ExperimentManager manager = new ExperimentManager(econf);
		
		SmithWatermanAlignment profile = new SmithWatermanAlignment(gconf, econf, manager); 
		
		ArgParser ap = new ArgParser(args);		
		int win = Args.parseInteger(args, "win", 200);
		
		List<Point> points = RegionFileUtilities.loadPeaksFromPeakFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		profile.setPoints(points);
		profile.setWidth(win);
		profile.loadData();
	
	}

}
