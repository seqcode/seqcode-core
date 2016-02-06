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
 * HyperstackSegmentation
 *
 * Methods for MultiScaleSignalRepresentation
 * Probabilistic Multiscale Image Segmentation, Vincken et al. IEEE (1997)
 * Multiscale representation of genomic signals, Knijnenburg et al. Nature Methods (2014)
 * 
 * @author naomi yamada
 * 
 * copied from SegmentationTree to take more than one experiment target
 *
 **/

public class HyperstackSegmentation {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected SEEDConfig sconfig;
	protected int numScale;
	
	/*********************
	 * Gaussian scale space and window parameters	
	 */
	final static double DELTA_TAU = 0.5*Math.log(2);
	final static double MINIMUM_VALUE = Math.pow(10, -100); //arbitrary minimum value; I cannot use Double.MIN_VALUE because it becomes zero after several multiplication with Double.MIN_VALUE
//	final static double P_MIN = Math.pow(10,-3); // Knijnenburg paper uses this value
	final static double P_MIN = Math.pow(10,-2); //needs to be carefully determined because it will substantially affect Gaussian window size
	final static double P_BIN = 0.4995; // to determine window size of growing bin
	final static double K_MIN = 1/Math.sqrt(1-Math.exp(-2*DELTA_TAU));	
	final static double K_N = Math.ceil(K_MIN);

	/*********************
	 * Linkage weights
	 */
	final static double WEIGHT_I = 1.00;
	final static double WEIGHT_G = 0.0000001;
	final static double WEIGHT_M = 1000;
		
	public HyperstackSegmentation(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, int scale){	
		gconfig = gcon;
		econfig = econ;
		sconfig = scon;
		numScale = scale; 
	}	

	protected Map <Integer, Set<Integer>> buildTree (float[][][] gaussianBlur, Map <Integer, Integer> linkageMap, float[] DImax, int trailingZero, int zeroEnd){
		
		int currchromBinSize = gaussianBlur.length;
		int numConditions = gaussianBlur[0][0].length;
		float fchromBinSize = currchromBinSize;
		
		//cumGaussianBlur stores cumulative signal counts; this to be used to get the mean signal intensity
		float[][][] cumGaussianBlur = new float [gaussianBlur.length][2][numConditions];
		for (int i = 0; i<currchromBinSize; i++){
			for (int j = 0; j<2; j++){
				for (int k = 0 ; k < numConditions ; k++)
					gaussianBlur[i][j][k] = 0;
			}
		}
		for (int k = 0 ; k < numConditions ; k++){
			cumGaussianBlur[0][1][k] = gaussianBlur[0][1][k];
			for (int i = 1; i < currchromBinSize; i++){
				cumGaussianBlur[i][1][k] = cumGaussianBlur[i-1][1][k] + gaussianBlur[i][1][k];
			}
		}
		
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
			NormalDistribution normDistribution = new NormalDistribution(0.00,sigma[n]);
			//take inverse CDF based on the normal distribution using probability
			double inverseCDF = normDistribution.inverseCumulativeProbability(P_MIN);	
			double binInverseCDF = normDistribution.inverseCumulativeProbability(P_BIN);			
			int windowSize = (int) (-Math.round(inverseCDF)*2+1);	
			float binWindowSize = -Math.round(binInverseCDF)*2+1;
			
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
			
			// smoothing signals for numConditions times
			for (int k = 0 ; k < numConditions; k++){
				
				double polyCoeffi[] = new double [(int) Math.ceil(fchromBinSize/binWindowSize)];
				System.out.println("binWindowSize is "+binWindowSize+"\tpolyCoeffi length is "+(int) Math.ceil(fchromBinSize/binWindowSize));
				//copy from column[1] to column[0];this procedure need to be repeated for each iteration of scale
				// copy from column[1] to array to store polynomial coefficient
				for (int i = 0 ; i<currchromBinSize; i++){
					gaussianBlur[i][0][k]=gaussianBlur[i][1][k];
					polyCoeffi[(int) Math.floor(((float) i)/binWindowSize)] += gaussianBlur[i][1][k]/binWindowSize;
				}
				
				for (int i = 0; i < Math.ceil(fchromBinSize/binWindowSize); i++){
					if (polyCoeffi[i] == 0)
						polyCoeffi[i] = MINIMUM_VALUE;
				}

				PolynomialFunction poly1 = new PolynomialFunction(polyCoeffi);
				PolynomialFunction poly2 = new PolynomialFunction(normalizedWindow);
				PolynomialFunction polyMultiplication=poly1.multiply(poly2);
				double coefficients[]= polyMultiplication.getCoefficients();
		
				//taking mid point of polynomial coefficients			
				int coeffiMid = (int) Math.floor(((float) coefficients.length)/ 2.0);
		
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
					if (polyCoeffi.length % 2 == 0 && coefficients.length % 2 == 1)
						gaussianBlur[i][1][k] = (float) coefficients[(int) (coeffiMid-Math.floor((fchromBinSize/2-i)/binWindowSize))+1];
					else
						gaussianBlur[i][1][k] = (float) coefficients[(int) (coeffiMid-Math.floor((fchromBinSize/2-i)/binWindowSize))];
				}			
			}
			
			// printing region is sensitive to bin size; current bin size is 100 bp
			for (int k = 0 ; k < numConditions ; k ++){
				System.out.println("condition "+k);
				for (int i = 0 ; i < gaussianBlur.length; i++){
					if (gaussianBlur[i][0][k] > 0){
						System.out.println("signal is "+gaussianBlur[i][0][k]+" index is "+i);
					}
				}
//				for (int i = 0; i< 1151698;i += 20000)
//					System.out.println(gaussianBlur[i][0][k]+" : "+gaussianBlur[i][1][k]);
			}
			
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

			final long linkageLoopStart = System.currentTimeMillis();
				
			/***************
			 * Linkage Loop
			 * 
			 * 1st iteration: only consider intensity difference between parent and kid and connect to the ones with the minimum differences
			 * From 2nd iteration: consider ground volume = # of nodes that parents are linked to the kids and mean signal intensity difference between parent and kid
			 * From 3rd iteration: increase the weight of the ground voluem by 1e-7
			 * Vincken paper states that after 3-4 iteration, there would be no significant difference	
			 */	
			
			TreeMap<Integer, Integer> GV_parents = new TreeMap<Integer,Integer>(linkageMap);
			double C_intensity = 0; 
			double C_ground = 0; 
			double C_meanIntensity = 0;
			double GV_max = 0;		
			double affectionScore = 0;
			Integer prevKid = 0;
			double M_kid = 0;
			double[] M_parent = new double [numConditions];
			//updating ground volume and iterating to encourage convergence
			for (int counter = 0; counter<5; counter++){
				
				if (counter != 0){
					for (Integer parent : GV_parents.keySet()){
						if ( GV_parents.get(parent) > GV_max)
							GV_max = GV_parents.get(parent);
					}				
				}	
				
				// look for parents within the windowSize
				for (Integer kid : linkageMap.keySet()){
					
					//per kid, we determine ground volume of parent and mean intensity of parent
					if (counter ==0 || GV_max == 0){
						C_ground = 0.00;
						C_meanIntensity = 0.00;
					}else{ // I changed weights of ground volume significantly from the previous version so I don't know if this will work in encouraging fewer parents
						C_ground = WEIGHT_G*counter*GV_parents.get(linkageMap.get(kid))/GV_max; 
						for (int k = 0 ; k < numConditions ; k++)
							// calculating the mean parents intensity; I need to double check this
							// line 252 is producing array out of bound exception
							M_parent[k] = cumGaussianBlur[linkageMap.get(kid)][0][k] - gaussianBlur[linkageMap.get(kid) - GV_parents.get(linkageMap.get(kid))][0][k];
					}
				
					double maxAffectionScore = 0;		

					for (int i = 0; i<DCPsize; i++){
						
						if (GV_parents.containsKey(kid+dcp[i])){
							
							double totalIdiff = 0;
							double totalMeanIdiff = 0;
							for (int k = 0 ; k <numConditions; k++){
								totalIdiff += Math.abs(gaussianBlur[kid][0][k] - gaussianBlur[kid+dcp[i]][1][k])/DImax[k];	
								M_kid = cumGaussianBlur[kid+dcp[i]][1][k] - cumGaussianBlur[prevKid][1][k];
								totalMeanIdiff += Math.abs(M_parent[k] - M_kid)/DImax[k];
							}
							C_intensity = WEIGHT_I*(1 - totalIdiff);
							C_meanIntensity = WEIGHT_M*(1-totalMeanIdiff);	

							//applying equation 9 in Vincken(1997)
							affectionScore = distanceFactor[i]*(C_intensity+C_ground+C_meanIntensity)/(WEIGHT_I + WEIGHT_G + WEIGHT_M);
							
							if (affectionScore > maxAffectionScore){
								maxAffectionScore = affectionScore;
								linkageMap.put(kid,(kid+dcp[i]));
							}							
						}							
					}
					prevKid = kid;
				}						
				//test
//						if (currchromBinSize > 20000000){			
//							System.out.println("current Chrom is: "+currChrom.getChrom());
//							System.out.println("printing linkangeMap content");
//							for (Map.Entry<Integer, Integer> entry : linkageMap.entrySet()){
//								System.out.println("Key: "+entry.getKey()+" Value: "+entry.getValue());
//							}
//						}
				GV_parents.clear();	
				Integer prevParent = 0;
				Map<Integer, Integer> sortedLinkageMap = new HashMap<Integer,Integer> (MapUtility.sortByValue(linkageMap));
				for (Integer parent : sortedLinkageMap.values()){
					GV_parents.put(parent, (parent-prevParent));
					prevParent = parent;
				}
				GV_parents.put(0, trailingZero);
				GV_parents.put(gaussianBlur.length-1,gaussianBlur.length-zeroEnd-1);
				
			}
			Map<Integer, Integer> sortedLinkageMap = new HashMap<Integer,Integer> (MapUtility.sortByValue(linkageMap));
			
			linkageMap.clear();
			for (Integer parent : sortedLinkageMap.values()){
				linkageMap.put(parent, parent);
			}						
			//for each scaleNum, add the parents to the segmentationTree
			
			final long linkageLoopEnd = System.currentTimeMillis();
			System.out.println("linkage Loop excusion time "+( linkageLoopEnd -linkageLoopStart));

			segmentationTree.put(n, GV_parents.keySet());
			
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

