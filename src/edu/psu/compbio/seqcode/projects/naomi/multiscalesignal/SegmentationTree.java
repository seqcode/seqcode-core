package edu.psu.compbio.seqcode.projects.naomi.multiscalesignal;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.distribution.NormalDistribution;

import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.projects.naomi.utilities.MapUtility;
import edu.psu.compbio.seqcode.projects.seed.SEEDConfig;

/**
 * Segmentation Tree
 *
 * Methods for MultiScaleSignalRepresentation
 * Probabilistic Multiscale Image Segmentation, Vincken et al. IEEE (1997)
 * 
 * @author naomi yamada
 *
 **/

public class SegmentationTree {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected SEEDConfig sconfig;
	protected int numScale;
	
	/*********************
	 * Gaussian scale space and window parameters	
	 */
	final static double DELTA_TAU = 0.5*Math.log(2);
	final static double MINIMUM_VALUE = Math.pow(10, -100); //arbitrary minimum value; I cannot use Double.MIN_VALUE because it can become zero
	// I have to determine P_MIN value carefully because P_MIN will substantially affect Gaussian window size
	final static double P_MIN = Math.pow(10,-3);
	//P_BIN is to determine window size of growing bin
	final static double P_BIN = 0.499;
	final static double K_MIN = 1/Math.sqrt(1-Math.exp(-2*DELTA_TAU));	
	final static double K_N = Math.ceil(K_MIN);

	/*********************
	 * Linkage parameters
	 */
	final static double WEIGHT_I = 1.00;
	final static double WEIGHT_G = 0.0000001;
	final static double WEIGHT_M = 1000;
		
	public SegmentationTree(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, int scale){	
		gconfig = gcon;
		econfig = econ;
		sconfig = scon;
		numScale = scale; 
	}	

	protected Map <Integer, Set<Integer>> buildTree (int currchromBinSize, float[][] gaussianBlur, Map <Integer, Integer> linkageMap, float DImax, int trailingZero, int zeroEnd){
		
		Map<Integer,Set<Integer>> segmentationTree =new HashMap<Integer,Set<Integer>>();
		segmentationTree.put(0, linkageMap.keySet());
		System.out.println("curr Scale 0 size, printing from segmentationTree "+segmentationTree.get(0).size());
		Set<Integer> startingNodes = new TreeSet<Integer>(segmentationTree.get(0));
		
//		Map<Region, HashMap<Integer,Set<Integer>>> segmentationTree = new HashMap<Region, HashMap<Integer, Set<Integer>>>();
		
		/*********************
		 * Matrices parameters
		 */	
		double[] sigma = new double[numScale];
		double[] radius = new double[numScale];
		for (int i = 0; i<numScale;i++){
			sigma[i] = 1;
			radius[i] = 1;
		}
		
		for (int n = 1; n<numScale; n++){
			
			final long gaussianStartTime = System.currentTimeMillis();
			
			//sigma calculation
			sigma[n] = Math.exp(n*DELTA_TAU);
			// create normal distribution with mean zero and sigma[n]
			NormalDistribution normDistribution = new NormalDistribution(0.00,sigma[n]);
			//take inverse CDF based on the normal distribution using probability
			double inverseCDF = normDistribution.inverseCumulativeProbability(P_MIN);	
			double binInverseCDF = normDistribution.inverseCumulativeProbability(P_BIN);			
			int windowSize = (int) (-Math.round(inverseCDF)*2+1);	
			int binWindowSize = (int) (-Math.round(binInverseCDF)*2+1);
			
			//window calculation based on Gaussian(normal) distribution with sigma, mean=zero,x=X[i]			
			double window[] = new double[windowSize];
			double windowSum = 0;
			for (int i = 0;i<windowSize;i++){
				window[i] = normDistribution.density(Math.round(inverseCDF)+i);
				windowSum = windowSum+window[i];
			}
			double normalizedWindow[]=new double[windowSize];
			for (int i = 0;i<windowSize;i++)
				normalizedWindow[i] = window[i]/windowSum;	
					
			double polyCoeffi[] = new double [(int) Math.ceil(currchromBinSize/binWindowSize)];
			System.out.println("binWindowSize is "+binWindowSize+"\tpolyCoeffi length is "+(int) Math.ceil(currchromBinSize/binWindowSize));
			//copy from column[1] to column[0];this procedure need to be repeated for each iteration of scale
			// copy from column[1] to array to store polynomial coefficient
			for (int i = 0 ; i<currchromBinSize; i++){
				gaussianBlur[i][0]=gaussianBlur[i][1];
				polyCoeffi[(int) Math.floor(i/binWindowSize)] += (gaussianBlur[i][1])/binWindowSize;
			}
			for (int i = 0; i < Math.ceil(currchromBinSize/binWindowSize); i++){
				if (polyCoeffi[i] == 0)
					polyCoeffi[i]=MINIMUM_VALUE;
			}	

			PolynomialFunction poly1 = new PolynomialFunction(polyCoeffi);
			PolynomialFunction poly2 = new PolynomialFunction(normalizedWindow);
			PolynomialFunction polyMultiplication=poly1.multiply(poly2);
			double coefficients[]= polyMultiplication.getCoefficients();
		
			//taking mid point of polynomial coefficients			
			int polyMid = (int) Math.floor(coefficients.length/2);
		
			System.out.println("currchromBin Size is : "+currchromBinSize+"\twindowSize is: "+windowSize+
					"\tpolyCoeffi length is "+polyCoeffi.length+"\tcoefficients length is: "+coefficients.length);

			//copy Gaussian blur results to the column[1] without increasing bin size
//			for (int i = 0; i<currchromBinSize;i++){
//				if (currchromBinSize % 2 ==0 && coefficients.length % 2 == 1){
//					gaussianBlur[i][1]=(float) coefficients[polyMid-currchromBinSize/2+i+1];
//				}else{
//					gaussianBlur[i][1]=(float) coefficients[polyMid-currchromBinSize/2+i];
//				}
			
			// copy Gaussian blur results to the column[1] with increasing bin size
			for (int i = 0; i<currchromBinSize;i++){
				gaussianBlur[i][1]=(float) coefficients[(int) (polyMid-Math.floor((currchromBinSize/2-i)/binWindowSize))];
				if (currchromBinSize % 2 ==0 && coefficients.length % 2 == 1)
					gaussianBlur[i][1]=(float) coefficients[(int) (polyMid-Math.floor((currchromBinSize/2-i-1)/binWindowSize))];
				else
					gaussianBlur[i][1]=(float) coefficients[(int) (polyMid-Math.floor((currchromBinSize/2-i)/binWindowSize))];
			}	
		
			for (int i = 0; i< 50;i++)
				System.out.println(gaussianBlur[11113388+i][0]+" : "+gaussianBlur[11113388+i][1]);
			
			final long gaussianEndTime = System.currentTimeMillis();
			
			System.out.println("scale "+n+" Gausisian blur execusion time "+ (gaussianEndTime-gaussianStartTime));
		
			/***************
			 * Search Volume
			 */ 	
		 
			double tempRadius;
			if (n==1){
				tempRadius = sigma[n];
			}else{
				tempRadius = Math.sqrt(Math.pow(sigma[n],2)-Math.pow(sigma[n-1], 2));
			}
			radius[n] = Math.ceil(K_MIN*tempRadius);
			
			int DCPsize = (int) (Math.round(radius[n])*2+1);
			int dcp[] = new int[DCPsize];
			double distanceFactor[] = new double[DCPsize];
			double affectionDistance;
			double denom = -2*(Math.pow(sigma[n], 2)-Math.pow(sigma[n-1],2));
		
			for (int i = 0; i<DCPsize;i++){
				dcp[i] = (int) -Math.round(radius[n])+i;
				// applying equation 7 in Vincken(1997)
				affectionDistance=Math.exp(Math.pow(dcp[i], 2)/denom)/Math.exp(Math.pow(0.5*sigma[n],2)/denom);
			
				//applying equation 8 in Vincken (1997) 
				if (Math.abs(dcp[i]) > 0.5*sigma[n]){distanceFactor[i]= affectionDistance;}
				else{distanceFactor[i] = 1.0000;}
			}
		
			/***************
			 * Linkage Loop	
			 */		 
			
			final long linkageLoopStart = System.currentTimeMillis();
			
//			TreeMap<Integer, Integer> GvParents = new TreeMap<Integer,Integer>();		
			TreeMap<Integer, Integer> GvParents = new TreeMap<Integer,Integer>(linkageMap);				 
			//First iteration only consider intensity differences between parent and kid and connect to the ones with the least difference.
			//From the second iteration, we consider ground volume = number of nodes that parents are linked to the kids
			//From third iteration, we increase the weight of the ground volume by 1e-7.
			//Vincken paper said after 3-4 iteration, there would be no significant difference.
			double groundVC = 0; 
			double groundVPmax = 0;		
			double tempScore = 0;
			//updating ground volume and iterating to encourage convergence
			for (int counter = 0; counter<5; counter++){
				
				if (counter != 0){
					for (Integer parent : GvParents.keySet()){
						if ( GvParents.get(parent) > groundVPmax)
							groundVPmax = GvParents.get(parent);
					}				
				}	
				
				// look for parents within the windowSize
				for (Integer kid : linkageMap.keySet()){
					
					if (counter ==0 || groundVPmax == 0){groundVC = 0.00;}
					else{ groundVC = (WEIGHT_I+WEIGHT_G*counter)*GvParents.get(linkageMap.get(kid))/groundVPmax;}
				
					double intensityDiffScore = 0;		

					for (int i = 0; i<DCPsize; i++){
						
						if (GvParents.containsKey(kid+dcp[i])){

							tempScore = distanceFactor[i]*((1- Math.abs(gaussianBlur[kid][0] - gaussianBlur[kid+dcp[i]][1])/DImax)+groundVC);
							
							if (tempScore > intensityDiffScore){
								intensityDiffScore = tempScore;
								linkageMap.put(kid,(kid+dcp[i]));
//test								if (counter ==0){linkageMap.put(kid,(kid+dcp[i]));}
//test								else{
				//					if(GvParents.containsKey(kid+dcp[i])){linkageMap.put(kid,(kid+dcp[i]));}
//									if(linkageMap.containsValue(kid+dcp[i])){linkageMap.put(kid,(kid+dcp[i]));}
//test								}
							}							
						}							
					}
				}						
				//test
				//		if (currchromBinSize > 20000000){			
				//			System.out.println("current Chrom is: "+currChrom.getChrom());
				//			System.out.println("printing linkangeMap content");
				//			for (Map.Entry<Integer, Integer> entry : linkageMap.entrySet()){
				//				System.out.println("Key: "+entry.getKey()+" Value: "+entry.getValue());
				//			}
				//		}
				GvParents.clear();	
				Integer lastParent = 0;
				Map<Integer, Integer> sortedLinkageMap = new HashMap<Integer,Integer> (MapUtility.sortByValue(linkageMap));
				for (Integer parent : sortedLinkageMap.values()){
					GvParents.put(parent, (parent-lastParent));
					lastParent = parent;
				}
				GvParents.put(0, trailingZero);
				GvParents.put(gaussianBlur.length-1,gaussianBlur.length-zeroEnd-1);
				
			}
			Map<Integer, Integer> sortedLinkageMap = new HashMap<Integer,Integer> (MapUtility.sortByValue(linkageMap));
			
			linkageMap.clear();
			for (Integer parent : sortedLinkageMap.values()){
				linkageMap.put(parent, parent);
			}						
			//for each scaleNum, add the parents to the segmentationTree
			
			final long linkageLoopEnd = System.currentTimeMillis();
			System.out.println("linkage Loop excusion time "+( linkageLoopEnd -linkageLoopStart));

			segmentationTree.put(n, GvParents.keySet());
			
		}//end of scale space iteration
		
		// scale zero is getting overwriting with the parents of the last scale; I'm overwriting the scale zero with initial nodesest for quick fix
		segmentationTree.put(0, startingNodes);
		
//		for (Integer scale : segmentationTree.keySet()){
//			System.out.println("current scale is: "+scale);
//			Set<Integer> sortedNodeSet = new TreeSet<Integer>(segmentationTree.get(scale));
//			System.out.println("current nodeset size is: "+sortedNodeSet.size());
//			for (Integer node : sortedNodeSet)
//				System.out.println(node);
//		}	
		
		return segmentationTree;
		
	}
}

