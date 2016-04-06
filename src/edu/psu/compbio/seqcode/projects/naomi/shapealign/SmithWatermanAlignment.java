package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

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
 * SmithWatermanAlignment: Aligns ChIP-exo reads using affine gaped Smith Waterman local alignment
 * 
 * @author naomi yamada
 *
 */

public class SmithWatermanAlignment {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	
	protected List<Point> points;
	protected List<Region> regions;
	protected int window;
	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;
	
	static final double MINIMUM_VALUE = -10000;
	
	protected Map<Sample, Map<Region,double[][]>> countsArray = new HashMap<Sample,Map<Region,double[][]>>();
	
	public SmithWatermanAlignment(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
	}
	
	// setters
	public void setPoints(List<Point> p){points = p;}
	public void setRegions(List<Region> reg){regions = reg;} 
	public void setWidth(int w){window = w;}
	public void setCountsArray(Map<Sample, Map<Region,double[][]>> sampleCounts){countsArray = sampleCounts;}
	
	public double computeScore(double aVal, double bVal){
		double score = (aVal + bVal)/2 - Math.abs(aVal - bVal);
		return score;		
	}
	
	public void loadData(){
		
		List<Region> region = new ArrayList<Region>();
		for(Point p: points){
			region.add(p.expand(window/2));
		}
		
		//merge overlapping regions and set back to the original size				
		List<Region> resizedRegions = new ArrayList<Region>();		
		for (Region reg : Region.mergeRegions(region)){ //mergeRegions is a static method	
			resizedRegions.add(reg.resize(window-1));
		}		
		setRegions(resizedRegions);

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
				double[][] sampleCounts = new double[window][2];
				for (int i = 0;i < window;i++){
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
	
	public void excuteShapeAlign(){
		
		for (Sample sample : countsArray.keySet()){
//			for (int i = 0; i <regions.size();i++){		
			
//			for (int i = 0; i <1 ; i++){	
			for (int i = 1; i <2 ; i++){
				for (int j = i+1; j <regions.size();j++)
					smithWatermanAlgorithm(sample, regions.get(i), regions.get(j));				
			}
		}				
	}
	
	public void smithWatermanAlgorithm(Sample sample, Region regA, Region regB){
		
		//get counts
		double [][] regACounts = countsArray.get(sample).get(regA);
		double [][] regBCounts = countsArray.get(sample).get(regB);
		
		//normalize the arrays to set the max value 1
		//should I be normalizing using max of either strand ?
		double maxA = MINIMUM_VALUE;
		double maxB = MINIMUM_VALUE;
		for (int i = 0; i <window ; i++){
			for (int s = 0 ; s < 2 ; s++){
				if (regACounts[i][s] > maxA){maxA = regACounts[i][s];}
				if (regBCounts[i][s] > maxB){maxB = regBCounts[i][s];}
			}
		}
		
		System.out.println("max counts are "+maxA+" : "+maxB);
		
		double [][] normRegACounts = new double [window][2];
		double [][] normRegBCounts = new double [window][2];
		double [][] normRegBRevCounts = new double[window][2];
		for (int i = 0 ; i <window; i++){
			for (int s = 0 ; s < 2 ; s++){
				normRegACounts[i][s] = regACounts[i][s]/maxA;
				normRegBCounts[i][s] = regBCounts[i][s]/maxB;
			}
		}
		
		//reversing normRegBCounts
		for (int i = 0; i <window; i++){
			for (int s = 0; s < 2 ;s++){
				normRegBRevCounts[window-i-1][1-s] = normRegBCounts[i][s];
			}
		}

		// align using two possible ways
		SmithWatermanAlignmentMatrix alignOne = new SmithWatermanAlignmentMatrix(normRegACounts,normRegBCounts);
		SmithWatermanAlignmentMatrix alignTwo = new SmithWatermanAlignmentMatrix(normRegACounts,normRegBRevCounts);
		alignOne.buildMatrix();
		alignTwo.buildMatrix();
		
		Stack<Integer> traceBack = new Stack<Integer>();
		double[][] regBarray = new double[window][2];
		boolean reverseB;
		int s_x_coord = 0;
		int s_y_coord = 0;
		int e_x_coord = 0;
		int e_y_coord = 0;
		
		if (alignOne.getMaxScore() > alignTwo.getMaxScore()){	
			
			regBarray = normRegBCounts;
			traceBack = alignOne.getTraceBack();
			s_x_coord = alignOne.getStartX();
			s_y_coord = alignOne.getStartY();
			e_x_coord = alignOne.getEndX();
			e_y_coord = alignOne.getEndY();		
			reverseB = false;
			
		}else{	
			
			regBarray = normRegBRevCounts;
			traceBack = alignTwo.getTraceBack();
			s_x_coord = alignTwo.getStartX();
			s_y_coord = alignTwo.getStartY();
			e_x_coord = alignTwo.getEndX();
			e_y_coord = alignTwo.getEndY();	
			reverseB = true;
		}
		
		System.out.println("alignment start coordinates "+ s_x_coord + " : " + s_y_coord);
		System.out.println("alignment end coordinates "+ e_x_coord + " : " + e_y_coord);
		
		System.out.println(Arrays.toString(traceBack.toArray()));
		
		double[][] alignedRegA = new double[e_x_coord-s_x_coord+1][2];
		double[][] alignedRegB = new double[e_y_coord-s_y_coord+1][2];
		int current_x = e_x_coord-1;
		int current_y = e_y_coord-1;
		
		// need to re-visit here. this is not working
		
		Stack<Integer> xTraceBack = traceBack;
		
		for (int i = e_x_coord-s_x_coord; i >= 0 ; i--){	
			
			for (int s = 0 ; s <2; s++)
				alignedRegA[i][s] = normRegACounts[current_x][s];
			
			if ( !xTraceBack.empty() ){
				if (xTraceBack.peek() == DIAG || xTraceBack.peek() == LEFT){
				current_x --;
				xTraceBack.pop();
				}
			}				
		}
		
		
		Stack<Integer> yTraceBack = traceBack;
		
		for (int i = e_y_coord-s_y_coord; i >= 0 ; i--){	
			
			for (int s = 0 ; s <2; s++){
				if (reverseB == false){
					alignedRegB[i][s] = normRegBCounts[current_y][s];
				}else{
					alignedRegB[i][s] = normRegBRevCounts[current_y][s];
				}
			}
			if ( !yTraceBack.empty() ){
				if (yTraceBack.peek() == DIAG || yTraceBack.peek() == UP){
				current_y --;	
				yTraceBack.pop();
				}
			}
		}
		
		System.out.println("before alignment reg A");
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(normRegACounts[i][0]+",");
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(normRegACounts[i][1]+",");	
		System.out.println();
		
		System.out.println("before alignment reg B");
		for (int i = 0; i < regBarray.length; i++)
			System.out.print(regBarray[i][0]+",");
		for (int i = 0; i < regBarray.length; i++)
			System.out.print(regBarray[i][1]+",");	
		System.out.println();
		
		System.out.println("aligned reg A");
		for (int i = 0; i < alignedRegA.length; i++)
			System.out.print(alignedRegA[i][0]+",");
		for (int i = 0; i < alignedRegA.length; i++)
			System.out.print(alignedRegA[i][1]+",");	
		System.out.println();

		System.out.println("aligned reg B");
		for (int i = 0; i < alignedRegB.length; i++)
			System.out.print(alignedRegB[i][0]+",");
		for (int i = 0; i < alignedRegB.length; i++)
			System.out.print(alignedRegB[i][1]+",");	
		System.out.println();
		
	}
	
	
	public static void main(String[] args){
		
		/***
		* Must be run with --fixedpb option; set it to high number >1000
		*
		***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		ExperimentManager manager = new ExperimentManager(econf);
			
//		econf.setPerBaseReadLimit(100000);
		econf.setPerBaseReadFiltering(false);
		
		SmithWatermanAlignment profile = new SmithWatermanAlignment(gconf, econf, manager); 
		
		ArgParser ap = new ArgParser(args);		
		int win = Args.parseInteger(args, "win", 200);
		
		List<Point> points = RegionFileUtilities.loadPeaksFromPeakFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		profile.setPoints(points);
		profile.setWidth(win);
		profile.loadData();
		profile.excuteShapeAlign();
	
	}

}
