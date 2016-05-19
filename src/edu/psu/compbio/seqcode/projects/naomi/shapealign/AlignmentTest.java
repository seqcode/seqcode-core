package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.projects.naomi.FeatureCountsLoader;

/**
 * Alignment Test: derived from SmithWaterman Alignment.java 
 * This code does not merge gff for near by regions
 * 
 * @author naomi yamada
 *
 */

public class AlignmentTest {
	
	protected FeatureCountsLoader featureCountsLoader;
	protected SimilarityScore similarity_s;
		
	protected List<StrandedRegion> strandeRegions;
	protected int window;
	
	protected double error = 0;
	protected double totalNum = 0;
	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;
	
	static final double MINIMUM_VALUE = -10000;
	
	protected Map<ControlledExperiment, Map<StrandedRegion, double[][]>> strandedRegionSampleCounts = new HashMap<ControlledExperiment, Map<StrandedRegion,double[][]>>();
	protected double [] offsetArray;
	
	public AlignmentTest(FeatureCountsLoader fcLoader, SimilarityScore sc){	
		featureCountsLoader = fcLoader;
		similarity_s = sc;
	}
	
	// setters
	public void setWidth(int w){window = w;}	
	
	// prints error rate
	public void printErrorRate(){System.out.println("error is "+error+ " totalNum is "+ totalNum+ " error rate is "+ (error/totalNum));}
	
	public void excuteShapeAlign(){
		
		strandedRegionSampleCounts = featureCountsLoader.strandedRegionSampleCounts();
		strandeRegions = featureCountsLoader.getStrandedRegions();
		
		//initiazlie offsetArray
		offsetArray = new double [window+1];
		for (int i = 0 ; i <= window ; i++)
			offsetArray[i] = 0;
		
		for (ControlledExperiment cExpt : strandedRegionSampleCounts.keySet()){
			for (int i = 0; i <strandeRegions.size();i++){		
			
//			for (int i = 0; i <1 ; i++){	
				for (int j = i+1; j <strandeRegions.size();j++){
//					System.out.println("region is "+regions.get(j).getLocationString());
					smithWatermanAlgorithm(cExpt, strandeRegions.get(i), strandeRegions.get(j));		
				}
			}
		}				
	}
	
	public void smithWatermanAlgorithm(ControlledExperiment cExpt, Region regA, Region regB){
		
		//get midpoints
//		double regAmid = regA.getMidpoint().getLocation();
//		double regBmid = regB.getMidpoint().getLocation();
		
		//get counts
		double [][] regACounts = strandedRegionSampleCounts.get(cExpt).get(regA);
		double [][] regBCounts = strandedRegionSampleCounts.get(cExpt).get(regB);
		
		//normalize the arrays to set the max value 1
		double maxA = MINIMUM_VALUE;
		double maxB = MINIMUM_VALUE;
		for (int i = 0; i <=window ; i++){
			for (int s = 0 ; s < 2 ; s++){
				if (regACounts[i][s] > maxA){maxA = regACounts[i][s];}
				if (regBCounts[i][s] > maxB){maxB = regBCounts[i][s];}
			}
		}
		
//		System.out.println("max counts are "+maxA+" : "+maxB);
		
		double [][] normRegACounts = new double [window+1][2];
		double [][] normRegBCounts = new double [window+1][2];
		double [][] normRegBRevCounts = new double[window+1][2];
		for (int i = 0 ; i <= window; i++){
			for (int s = 0 ; s < 2 ; s++){
				normRegACounts[i][s] = regACounts[i][s]/maxA;
				normRegBCounts[i][s] = regBCounts[i][s]/maxB;
			}
		}
		
		//reversing normRegBCounts
		for (int i = 0; i <= window; i++){
			for (int s = 0; s < 2 ;s++){
				normRegBRevCounts[window-i][1-s] = normRegBCounts[i][s];
			}
		}

		// align using two possible ways
		SmithWatermanAlignmentMatrix alignOne = new SmithWatermanAlignmentMatrix(normRegACounts,normRegBCounts, similarity_s);
		SmithWatermanAlignmentMatrix alignTwo = new SmithWatermanAlignmentMatrix(normRegACounts,normRegBRevCounts, similarity_s);
		alignOne.buildMatrix();
		alignTwo.buildMatrix();
		
		Stack<Integer> traceBack = new Stack<Integer>();
		double[][] regBarray = new double[window+1][2];
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
		
		double x_mid = (s_x_coord + e_x_coord)/2;
		double y_mid = (s_y_coord + e_y_coord)/2;
		
//		System.out.println("alignment start coordinates "+ s_x_coord + " : " + s_y_coord);
//		System.out.println("alignment end coordinates "+ e_x_coord + " : " + e_y_coord);
		
//		System.out.println(Arrays.toString(traceBack.toArray()));
		
		double[][] alignedRegA = new double[e_x_coord-s_x_coord+1][2];
		double[][] alignedRegB = new double[e_y_coord-s_y_coord+1][2];
		int current_x = e_x_coord-1;
		int current_y = e_y_coord-1;
		
		// trace back for x 
		@SuppressWarnings("unchecked")
		Stack<Integer> xTraceBack = (Stack<Integer>) traceBack.clone();
		
		for (int i = e_x_coord-s_x_coord; i >= 0 ; i--){	
			
			if (current_x >= 0){

				for (int s = 0 ; s <2; s++)
					alignedRegA[i][s] = normRegACounts[current_x][s];
			
				if ( !xTraceBack.empty() ){			
					if (xTraceBack.peek() == DIAG || xTraceBack.peek() == LEFT){
						current_x --;
					}			
					xTraceBack.pop();
				}	
			}
		}		
		
		// trace back for y 
		@SuppressWarnings("unchecked")
		Stack<Integer> yTraceBack = (Stack<Integer>) traceBack.clone();;
		
		for (int i = e_y_coord-s_y_coord; i >= 0 ; i--){	
			
			if (current_y >= 0){
			
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
					}
					yTraceBack.pop();			
				}
			}
		}	

		/**
		System.out.println("before alignment reg A");
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(normRegACounts[i][0]+",");
		System.out.println();
		for (int i = 0; i < normRegACounts.length; i++)
			System.out.print(normRegACounts[i][1]+",");	
		System.out.println();
		
		System.out.println("before alignment reg B");
		for (int i = 0; i < normRegBCounts.length; i++)
			System.out.print(normRegBCounts[i][0]+",");
		System.out.println();
		for (int i = 0; i < normRegBCounts.length; i++)
			System.out.print(normRegBCounts[i][1]+",");	
		System.out.println();
		
		System.out.println("aligned reg A");
		for (int i = 0; i < alignedRegA.length; i++)
			System.out.print(alignedRegA[i][0]+",");
		System.out.println();
		for (int i = 0; i < alignedRegA.length; i++)
			System.out.print(alignedRegA[i][1]+",");	
		System.out.println();

		System.out.println("aligned reg B");
		for (int i = 0; i < alignedRegB.length; i++)
			System.out.print(alignedRegB[i][0]+",");
		System.out.println();
		for (int i = 0; i < alignedRegB.length; i++)
			System.out.print(alignedRegB[i][1]+",");	
			
		**/
		
		// increment offset array
		offsetArray[(int) (y_mid-x_mid+window/2)]++;				
		// incrementing error allowing offset of +-1
		totalNum += 1;
		if ( traceBack.contains(LEFT) || traceBack.contains(UP) ){ // check that stack only contains DIAG
			error += 1;
		}else{
			if (Math.abs(y_mid-x_mid) >1){
				error +=1;
			}
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
		
		FeatureCountsLoader fcLoader = new FeatureCountsLoader(gconf, econf, manager);
		fcLoader.setStrandedPoints(spoints);
		fcLoader.setWindowSize(win); // window size must be twice bigger so it can slide window size times

		SimilarityScore sc = new SimilarityScore();
		
		if (Args.parseFlags(args).contains("euclidean")){ sc.setEuclideanL2();}
		if (Args.parseFlags(args).contains("sorensen")){ sc.setSorensen();}
		if (Args.parseFlags(args).contains("soergel")){ sc.setSoergel();}
		if (Args.parseFlags(args).contains("lorentzian")){ sc.setLorentzian();}
		if (Args.parseFlags(args).contains("pce")){ sc.setPCE();}
		if (Args.parseFlags(args).contains("squaredChi")){ sc.setSquaredChi();}
		if (Args.parseFlags(args).contains("divergence")){ sc.setDivergence();}
		if (Args.parseFlags(args).contains("clark")){ sc.setClark();}
		

		AlignmentTest profile = new AlignmentTest(fcLoader, sc); 		
		profile.setWidth(win);
		profile.excuteShapeAlign();
		profile.printErrorRate();	
		if (Args.parseFlags(args).contains("printOffset")){profile.printOffsetArray();}
		
		manager.close();
	}
}
