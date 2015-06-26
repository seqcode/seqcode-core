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
//import edu.psu.compbio.seqcode.gse.utils.ArgParser;

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
	
	public void K_LineMeans(int K){
		
		System.out.println("current K is: "+K);
		
		List<Double> angles = new ArrayList<Double>();	
		//fix when k==1 . currently the program cannot handle when k==1
		if (K==1){
			angles.add((double) 45);
		}else if (K>=2){
			angles.add((double)0.000000000000001);
			double seg = 90/(K-1);
			System.out.println("seg value is: "+seg);
			if (seg<90){
				while (seg<90){
					angles.add(seg);
					seg+=seg;
				}
			}
			angles.add((double) 90-0.000000000000001);
		}
		
		//testing 
		System.out.println("contents of angle: ");
		for (int i = 0; i<angles.size();i++)
			System.out.println(angles.get(i));
		
		// converting angle in degree to radians, and tangent to get slopes
		List<Double> slopes = new ArrayList<Double>();
		for (double angle : angles)
			slopes.add(Math.tan(Math.toRadians(angle)));
		
		//testing
		System.out.println("contents and index of slopes: ");
		for (double slope: slopes)
			System.out.println(slope+"\t"+slopes.indexOf(slope));
		
		angles.clear();
		
		List<Double> previousSlopes = new ArrayList<Double>();
		for (int i = 0; i<slopes.size();i++)
			previousSlopes.add((double) 0);
		
		//testing
		System.out.println("contents of previousSlopes: ");
		for (int i = 0; i<previousSlopes.size();i++)
			System.out.println(previousSlopes.get(i));

		
		float [][] xyPairs = getXYpairs();		
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
		while (!previousSlopes.equals(slopes)||iteration_tracker<10){ 
			
			//testing
			System.out.println("**************inside while loop*******************");
			System.out.println("content of previousSlopes is: ");
			for (double previous : previousSlopes)
				System.out.println(previous);
			
			System.out.println("contents of slopes is: ");
			for (double sl : slopes)
				System.out.println(sl);
		
			//calculating intersect, intersecting x and y points and squared distances from slope
			double neg_inverse = 0;
			double intersect = 0;
			double x_point = 0;
			double y_point = 0;
			double squareDistance = 0;
			for (double slope : slopes){
				neg_inverse=-(1/slope);
				
				//testing
				System.out.println("neg_inverse is: "+neg_inverse);
				
				for (int i = 0; i<xyPairs.length; i++){
					intersect = xyPairs[i][1]-(neg_inverse*xyPairs[i][0]);
					x_point = intersect/(slope-neg_inverse);
					y_point = slope*x_point;
					squareDistance = Math.pow(x_point-xyPairs[i][0],2)+Math.pow(y_point-xyPairs[i][1],2);
					//error! array out of bound exception					
					distanceArray[i][slopes.indexOf(slope)]= squareDistance;
					
					//initialize everything
					intersect = 0;
					x_point = 0;
					y_point = 0;
					squareDistance = 0;
				}
				neg_inverse = 0;
			}
						
			//testing
			System.out.println("contents of distance Array: ");
			for (int i = 0; i<10;i++){
				for (int j = 0; j <K; j++)
					System.out.println(distanceArray[i][j]);
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
					System.out.println("min value and index are: "+minimum+"\t"+minIndex);
				}
				xySlopes[i] = slopes.get(minIndex);	
				minimum = Integer.MAX_VALUE;
				minIndex = 0;
			}
			
			//testing
			System.out.println("contents of xySlopes: " );
			for (int i = 0; i <10;i++)
				System.out.println(xySlopes);
		
			float xSum = 0;
			float ySum = 0;
			int N = 0;	
			double xMeans = xSum/N;
			double yMeans = ySum/N;
			//copy slopes in previousSlopes and remove contents of slopes
			previousSlopes.clear();
			previousSlopes.addAll(slopes);
			
			//testing
			System.out.println("contents of previousSlopes: ");
			for (double s : previousSlopes)
				System.out.println(s);
			
			slopes.clear();
			
			//comparing the values of slope not the object; this may cause issues when excecuting
			for (double slope : previousSlopes){
				//testing
				System.out.println("current slope is: "+slope);
				
				for (int i = 0; i<xyPairs.length;i++){
					if (xySlopes[i]==slope){
						xSum+=xyPairs[i][0];
						ySum+=xyPairs[i][1];
						N++;
					}
				}
				xMeans = xSum/N;
				yMeans = ySum/N;
				slopes.add(yMeans/xMeans);
				
				System.out.println("xSum is: "+xSum+" ySum is: "+ySum+" N is "+N+" xMeans is "+xMeans+" yMeans is "+yMeans);
				
				//initialize everything
				xSum = 0;
				ySum = 0;
				N = 0;
				xMeans = 0;
				yMeans = 0;

			}
			
			System.out.println("current iteration number is: "+iteration_tracker);
			System.out.println("printing current list of slopes");
			for (double slope : slopes)
				System.out.println(slope);
			
			iteration_tracker++;
			
		}//finish the loop once the values in slopes stop changing
	}
	
	public static void main(String[] args){
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		
		//later add option to take different k 
//		ArgParser ap = new ArgParser(args);
				
		CrossContaminationEstimator estimator = new CrossContaminationEstimator (gconf, econf);
//		estimator.getXYpairs();	
//		estimator.printXYpairs();
		// trying with various k
		for (int k = 2; k<=7;k++)
			estimator.K_LineMeans(k);
		
	}
}