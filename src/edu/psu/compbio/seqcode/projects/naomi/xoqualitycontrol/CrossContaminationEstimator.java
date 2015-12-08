package edu.psu.compbio.seqcode.projects.naomi.xoqualitycontrol;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class CrossContaminationEstimator {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected float[][] xyPairs;


	public final float CONST1 = 2000000;
	
	public CrossContaminationEstimator(GenomeConfig gcon, ExptConfig econ){	
		gconfig = gcon;
		econfig = econ;
	}
	
	public void loadXYpairs(){
		
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
//		int dataPointsSize = 0;		
//		for (int i = 0; i<(int) genome.getGenomeLength(); i++){
//			if (dataPoints[i][0]!=0)
//				dataPointsSize++;
//		}		
		//copy non-zero dataPoints to xyPoints
//		float [][] xyPairs = new float[dataPointsSize][3];		
//		for (int i = 0; i<(int) genome.getGenomeLength();i++){
//			if (dataPoints[i][0]!=0){
//				for (int s= 0; s<3;s++)
//					xyPairs[i][s]=dataPoints[i][s];
//			}
//		}
		
		//only copying datapoints which go over some upper limits

		int SumAllCounts = 0;
		for (Sample sample: manager.getSamples()){
			SumAllCounts+= sample.getHitCount();
		}		
		double upperLimit = 0;
		upperLimit = SumAllCounts/CONST1;
		float CONST=10;
		
		System.out.println("upperLimit is: "+upperLimit);
		
		int dataPointSize = 0;
		for (int i = 0; i<(int) genome.getGenomeLength(); i++){
			if (((dataPoints[i][0]+dataPoints[i][1])>upperLimit)&&(dataPoints[i][0]+dataPoints[i][1])>CONST)
				dataPointSize++;
		}		
		System.out.println("dataPointSize is: "+dataPointSize);		
		int xy_index = 0;		
		xyPairs = new float[dataPointSize][3];
		for (int i = 0; i<(int) genome.getGenomeLength();i++){
			if (((dataPoints[i][0]+dataPoints[i][1])>upperLimit)&&(dataPoints[i][0]+dataPoints[i][1])>CONST){
				for (int s = 0; s<3;s++){
					xyPairs[xy_index][s]= dataPoints[i][s];
				}
				xy_index++;
			}				
		}	
	}

	public void printXYpairs(String out) throws FileNotFoundException, UnsupportedEncodingException {
		
		PrintWriter writer = new PrintWriter(out,"UTF-8");
		//printing xyPairs
		writer.println("#max_tag_number\tsum_of_other_sample's_tags\tsampleID");
		for (int i = 0; i< xyPairs.length; i++)
			writer.println(xyPairs[i][0]+"\t"+xyPairs[i][1]+"\t"+xyPairs[i][2]);
		writer.close();
	}			
	
	public void K_LineMeans(int K){
		
		System.out.println("*******computing slopes in K_LinesMeans********current K is: "+ K);
		
		List<Double> angles = new ArrayList<Double>();	
		
		//initializing angles according to number of K
		if (K==1){
			angles.add((double) 45);
		}else if (K>=2){
			angles.add((double)0.000000000000001);
			double seg = (double)90/((double)(K-1));
			double total = seg;
			System.out.println("seg value is: "+seg);
			if (total<90){
				while (total<90){
					angles.add(total);
					total+=seg;
				}
			}
			angles.add((double) 90-0.000000000000001);
		}
		
		// converting angle in degree to radians, and tangent to get slopes
		List<Double> slopes = new ArrayList<Double>();
		for (double angle : angles)
			slopes.add(Math.tan(Math.toRadians(angle)));
		
		angles.clear();
		
		List<Double> previousSlopes = new ArrayList<Double>();
		for (int i = 0; i<slopes.size();i++)
			previousSlopes.add((double) 0);
			
		//xyPairs_slope is a parallel double arrays of xyPairs (hence, same index system) that hold slope values.
		double [] xySlopes = new double [xyPairs.length];
		
		//testing
		System.out.println("xyPairs.length is: "+xyPairs.length);
		
		//distanceArray holds squared distance of each point to each slope
		double distanceArray[][] = new double [xyPairs.length][K];		
		
		//this is to test how many iteration it is making
		int iteration_tracker = 1;
		
		//iterating till slopes stop changing
		//error !iteration is not stopping by comparing to the previous slope list!!!!
		while (!previousSlopes.equals(slopes)){ 
			
			//calculating intersect, intersecting x and y points and squared distances from slope
			double neg_inverse = 0;
			double intersect = 0;
			double x_point = 0;
			double y_point = 0;
			double squareDistance = 0;
			for (double slope : slopes){
				neg_inverse=-(1/slope);
				
				for (int i = 0; i<xyPairs.length; i++){
					intersect = xyPairs[i][1]-(neg_inverse*xyPairs[i][0]);
					x_point = intersect/(slope-neg_inverse);
					y_point = slope*x_point;
					squareDistance = Math.pow(x_point-xyPairs[i][0],2)+Math.pow(y_point-xyPairs[i][1],2);
										
					distanceArray[i][slopes.indexOf(slope)]= squareDistance;
					
					//initialize everything back to zeros
					intersect = 0;
					x_point = 0;
					y_point = 0;
					squareDistance = 0;
				}
				neg_inverse = 0;
			}

			//find minimum distance and put the distance in xySlopes array
			double minimum = Integer.MAX_VALUE;
			int minIndex = 0;
			for (int i = 0; i <xyPairs.length;i++){
				for (int s = 0; s<K;s++){
					if (distanceArray[i][s]<minimum){
						minimum = distanceArray[i][s];
						minIndex = s;
					}
				}
				xySlopes[i] = slopes.get(minIndex);	
				minimum = Integer.MAX_VALUE;
				minIndex = 0;
			}
			
			//copy slopes in previousSlopes and remove contents of slopes
			previousSlopes.clear();
			previousSlopes.addAll(slopes);
			slopes.clear();
						
			//if the new slope is not equal to the old slope, re-calculate slopes
			//apply regression here
			for (double slope : previousSlopes){
				
				SimpleRegression regression = new SimpleRegression(false);
								
				for (int i = 0; i<xyPairs.length;i++){
					if (xySlopes[i]==slope){
						regression.addData(xyPairs[i][0],xyPairs[i][1]);
					}
				}
				slopes.add(regression.getSlope());
				
				System.out.println("number of observation is "+regression.getN()+" correlation is "+regression.getR()+" slope is "+regression.getSlope());
			}
			
			System.out.println("current iteration number is: "+iteration_tracker);
			System.out.println("printing current list of slopes");
			for (double slope : slopes)
				System.out.println(slope);
			
			iteration_tracker++;
			
		}//finish the loop once the values in slopes stop changing
	}
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException{
		
		/***
		 * You need to specify --fixedpb otherwise upper count limit would be set to a random number
		 * The method can take --count --file [file name] 
		 * or
		 * --K [number] or --varK
		 ***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		
		CrossContaminationEstimator estimator = new CrossContaminationEstimator (gconf, econf);
		
		estimator.loadXYpairs();
		
		ArgParser ap = new ArgParser(args);
		
		if (ap.hasKey("count") && ap.hasKey("file")){
			String outName = null;
			outName = Args.parseString(args, "file", "count.txt");
			estimator.printXYpairs(outName);
		}
		
		if (ap.hasKey("K")){
			int k_num = Args.parseInteger(args,"K",5);
			estimator.K_LineMeans(k_num);
		}else if (ap.hasKey("varK")){
			for (int k = 1;k<=4;k++)
				estimator.K_LineMeans(k);
		}			
	}

}
