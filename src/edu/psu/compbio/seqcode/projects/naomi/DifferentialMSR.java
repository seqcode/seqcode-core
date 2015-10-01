package edu.psu.compbio.seqcode.projects.naomi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.distribution.NormalDistribution;

import pal.util.Comparator;
import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.seed.SEEDConfig;

/**
 * DifferentialMSR: 
 *
 * Methods refer to two papers
 * Probabilistic Multiscale Image Segmentation, Vincken et al. IEEE (1997)
 * 
 * @author naomi yamada
 *
 **/

public class DifferentialMSR {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected SEEDConfig sconfig;
	
	protected int threePrimReadExt = 200;
	protected int binWidth = 1;
	
	//Map<chrm, Map<scale space, Set<segmentation>>>
	protected Map<Region, Map<Integer,Set<Integer>>> segmentationTree;

	public DifferentialMSR(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon){	
		gconfig = gcon;
		econfig = econ;
		sconfig = scon;
	}
	 
	public void buildMSR(){
		
		/*********************
		 * Gaussian scale space and window parameters	
		 */
		// arbitrary number of scale
		int numScale= 20;
		double DELTA_TAU = 0.5*Math.log(2);	
		double MINIMUM_VALUE = Math.pow(10, -100); //arbitrary minimum value; I cannot use Double.MIN_VALUE because it can become zero
		// I have to determine P_MIN value carefully because P_MIN will substantially affect Gaussian window size
		double P_MIN = Math.pow(10,-3);
		double K_MIN = 1/Math.sqrt(1-Math.exp(-2*DELTA_TAU));	
		double K_N = Math.ceil(K_MIN);
		
		/*********************
		 * Linkage parameters
		 */
		double WEIGHT_I = 1.00;
		double WEIGHT_G = 0.0000001;
		double WEIGHT_M = 1000;
		
		/*********************
		 * Matrices parameters
		 */	
		double sigma[] = new double[numScale];
		double radius[] = new double[numScale];
		for (int i = 0; i<numScale;i++){
			sigma[i] = 1;
			radius[i] = 1;
		}
		
		ExperimentManager manager = new ExperimentManager(econfig);
		Genome genome = gconfig.getGenome();
		
		//test to print whole chromosomes
		System.out.println(genome.getChromList());
		
		//fix here to get parameters only if they are specified
		binWidth = sconfig.getBinWidth();
		threePrimReadExt = sconfig.getTag3PrimeExtension();
		
		//test to print binWidth and threePrimReadExt
		System.out.println("binWidth is: "+binWidth);
		System.out.println("threePrimReadExt is: "+threePrimReadExt);
		
		Iterator<Region> chroms = new ChromosomeGenerator<Genome>().execute(genome);
		//iterating each chromosome (each chromosome is a region).
		while (chroms.hasNext()) {
			
			Region currChrom = chroms.next();	
			int currchromSize = currChrom.getWidth();
			int currchromBinSize = (int) Math.ceil(currchromSize/binWidth);
			
			//primitive matrix to store signal and the subsequent convolved signals
			//its index correspond to the coordinates
			float[][] GaussianBlur = new float[currchromBinSize][2];
			for (int i = 0;i<currchromBinSize;i++){
				for (int j = 0; j<2;j++)
					GaussianBlur[i][j] = 0;
			}
			
			//get StrandedBaseCount list for each chromosome
			Map<Sample, List<StrandedBaseCount>> sampleCountsMap = new HashMap<Sample, List<StrandedBaseCount>>();
			for (Sample sample : manager.getSamples())
				sampleCountsMap.put(sample,sample.getBases(currChrom)); 
			
			//StrandedBasedCount object contains positive and negative strand separately
			//store all base counts indexed by positions at column[1]
			//extend reads to 3' end and bin according to bin size			
			for (Sample sample : manager.getSamples()){		
				List<StrandedBaseCount> currentCounts = sampleCountsMap.get(sample);
				for (StrandedBaseCount hits: currentCounts){
					for (int i = 0; i<threePrimReadExt+1; i++){
						if (hits.getStrand()=='+' && hits.getCoordinate()+i<currchromSize){
								GaussianBlur[(int) Math.ceil((hits.getCoordinate()+i)/binWidth)][1]+=hits.getCount();
						}else if (hits.getStrand()=='+' && hits.getCoordinate()-i >=0){
								GaussianBlur[(int) Math.ceil((hits.getCoordinate()-i)/binWidth)][1]+=hits.getCount();
						}
					}
				}
				currentCounts = null;
			}
			
			//testing
//			if (currchromSize > 200000000){			
//				System.out.println("current Chrom is: "+currChrom.getChrom());
//				for (int i = 0; i< 100;i++)
//					System.out.println(GaussianBlur[(int) Math.ceil((92943501)/binWidth)+i][1]);
//			}
			
			/*********************
			 * Starting nodes
			 */
					
			//linkageMap contains index of kids and parents
			HashMap <Integer, Integer> linkageMap = new HashMap<Integer, Integer>();
			//adding starting nodes; to qualify for the starting nodes the signal intensity needs to be different from the subsequent signal intensity
			//adding the starting and end positions in the kids at start and end positions  
			//setting max & min signal intensity  
			float DImax = 0;
			float DImin = (float) Integer.MAX_VALUE;
			List <Integer> nonzeroList = new ArrayList<Integer>();
			linkageMap.put(0,0);
			for (int i = 0 ; i< GaussianBlur.length-1; i++){ //should I start
				if (GaussianBlur[i][1] != GaussianBlur[i+1][1])
					linkageMap.put(i,0);
				if (GaussianBlur[i][1] > DImax)
					DImax = GaussianBlur[i][1];
				if (GaussianBlur[i][1] < DImin)
					DImin = GaussianBlur[i][1];		
				if (GaussianBlur[i][1]!=0)
					nonzeroList.add(i);
			}
			linkageMap.put(GaussianBlur.length-1,0);
			
			//copy to segmentation tree
			
			Map<Integer,Set<Integer>> currScale =new HashMap<Integer,Set<Integer>>();
			currScale.put(1, linkageMap.keySet());
			System.out.println("curr Scale 1 size: "+linkageMap.keySet().size());
			
			//determine the first nonzero and last nonzero from signal
			
			//check here ; sometimes it produces zero for DImax, DImin, trailingZero, zeroEnd
			int trailingZero = 0;
			int zeroEnd = 0;		
			if (!nonzeroList.isEmpty()){
				trailingZero = Collections.min(nonzeroList)-1;
				zeroEnd = Collections.max(nonzeroList)+1;
			}
			if (trailingZero == -1)
				trailingZero = 0;
			
			System.out.println("DImax is: "+DImax+"\t"+"DImin is: "+DImin+
					"\t"+"trailingZero: "+trailingZero+"\t"+"zeroEnd"+"\t"+zeroEnd);		

			for (int n = 1;n < numScale; n++){
				
				/*********************
				 * Gaussian scale space 
				 */	
				 
				double polyCoeffi[] = new double [currchromBinSize];
				//first copy from column[1] to column[0];this procedure need to be repeated for each iteration of scale
				//also copy from column[1] to array to store polynomial coefficient
				for (int i = 0 ; i<currchromBinSize; i++){
					GaussianBlur[i][0]=GaussianBlur[i][1];
					if (GaussianBlur[i][1] != 0){
						polyCoeffi[i]=GaussianBlur[i][1];
					}else{
						polyCoeffi[i]=MINIMUM_VALUE;
					}
				}
				//sigma calculation
				sigma[n] = Math.exp(n*DELTA_TAU);
				// create normal distribution with mean zero and sigma[n]
				NormalDistribution normDistribution = new NormalDistribution(0.00,sigma[n]);
				//take inverse CDF based on the normal distribution using probability
				double inverseCDF = normDistribution.inverseCumulativeProbability(P_MIN);				
				int windowSize = (int) (-Math.round(inverseCDF)*2+1);						
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
	
				//multiplying by polynomial ; I have to test to see how this works
				PolynomialFunction poly1 = new PolynomialFunction(polyCoeffi);
				PolynomialFunction poly2 = new PolynomialFunction(normalizedWindow);
				PolynomialFunction polyMultiplication=poly1.multiply(poly2);
				double coefficients[]= polyMultiplication.getCoefficients();
				
				//taking mid point of polynomial coefficients			
				int polyMid = (int) Math.floor(coefficients.length/2);
				
				System.out.println("currchromBin Size is : "+currchromBinSize+"\t"+ "windowSize is: "+windowSize+"\t"+"coefficients length is: "+coefficients.length);
		
				//copy Gaussian blur results to the column[1]
				// I should check to make sure that it's not off by 1
				for (int i = 0; i<currchromBinSize;i++){
					if (currchromBinSize % 2 ==0 && coefficients.length/2 == 1)
						GaussianBlur[i][1]=(float) coefficients[polyMid-currchromBinSize/2+i+1];
					else
						GaussianBlur[i][1]=(float) coefficients[polyMid-currchromBinSize/2+i];
				}	
				
				//testing
//				if (currchromSize > 200000000){			
//					System.out.println("current Chrom is: "+currChrom.getChrom());
//					for (int i = 0; i< 100;i++)
//						System.out.println(GaussianBlur[(int) Math.ceil((92943501)/binWidth)+i][0]+" : "+GaussianBlur[(int) Math.ceil((92943501)/binWidth)+i][1]);
//				}
				
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
				 
				TreeMap<Integer, Integer> GvParents = new TreeMap<Integer,Integer>();			
					/***********
					 * Over window
					 */					 
				//build segmentTree 
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
							//307 is a problem because all parents are zero
							if (GvParents.containsValue(GvParents.get(parent))){
								if ( GvParents.get(parent) > groundVPmax)
									groundVPmax = GvParents.get(parent);
							}
							else {System.out.println("no value was found! ");}
						}				
					}	
						
					for (Integer key : linkageMap.keySet()){
						System.out.println("child is : "+key+" parents is: "+linkageMap.get(key));
					}
					for (Integer key : GvParents.keySet()){
						System.out.println("parents is : "+key);
					}
						
					System.out.println("groundVPmax is "+groundVPmax);
					for (Integer kid : linkageMap.keySet()){
						double intensityDiffScore = 0;
							
						//there is something wrong with iteration of DCPsize and updating parents
						//iteration is always off by 2. is it important to adjust for that?
							
						for (int i = 0; i<DCPsize; i++){
							if ((kid + dcp[i]) >=0 && (kid + dcp[i]) <currchromBinSize){
								if (counter ==0 || groundVPmax == 0){
									groundVC = 0.00;
								//do linkageMap.get(kid) and GvParents keys always match? i don't think so ; where are we modifying parents?
								// I can make linkageMap <kids<parent,groundVolume> or maybe parents are updated everytime?
								}else{ 		
//									Integer Gv;
//									if (GvParents.containsKey(linkageMap.get(kid))){Gv = GvParents.get(linkageMap.get(kid));}
//									else if(tempGvParents.containsKey(linkageMap.get(kid))){Gv = tempGvParents.get(linkageMap.get(kid));}
//									else {System.out.println("parents are missing !");
									groundVC = (WEIGHT_I+WEIGHT_G*counter)*GvParents.get(linkageMap.get(kid))/groundVPmax;
								}
									
								tempScore = distanceFactor[i]*((1- Math.abs(GaussianBlur[kid][0] - GaussianBlur[kid+dcp[i]][1])/DImax)+groundVC);
								if (tempScore > intensityDiffScore){
									intensityDiffScore = tempScore;
									GvParents.put((kid+dcp[i]),GvParents.get(linkageMap.get(kid))); //add parents in GvParents
										
//									GvParents.remove(linkageMap.get(kid)); //remove previous parents
									linkageMap.put(kid,(kid+dcp[i]));	//update parents in linkageMap	
//									System.out.println("intensityDiffScore is: "+intensityDiffScore+" DImax is "+DImax);
//									System.out.println("kid is: "+kid+" value is "+(kid+dcp[i]));
								}
								//updating GvParents
							}							
						}
					}
						
						//test
				//		if (currchromSize > 200000000){			
				//			System.out.println("current Chrom is: "+currChrom.getChrom());
				//			System.out.println("printing linkangeMap content");
				//			for (Map.Entry<Integer, Integer> entry : linkageMap.entrySet()){
				//				System.out.println("Key: "+entry.getKey()+" Value: "+entry.getValue());
				//			}
				//		}
												
					Integer lastParent = 0;
					Map<Integer, Integer> sortedLinkageMap = MapUtility.sortByValue(linkageMap);
					for (Integer parent : sortedLinkageMap.values()){
						GvParents.put(parent, (parent-lastParent));
						lastParent = parent;
					}
					GvParents.put(0, trailingZero);
					GvParents.put(GvParents.firstKey(), trailingZero-GvParents.firstKey());
	//					GvParents.put(GaussianBlur.length,GaussianBlur.length-zeroEnd);			
						
						//test
//						if (currchromSize > 200000000){			
//							System.out.println("current Chrom is: "+currChrom.getChrom());
//							System.out.println("printing GvParents content");
//							for (Map.Entry<Integer, Integer> entry : GvParents.entrySet()){
//								System.out.println("Key: "+entry.getKey()+" Value: "+entry.getValue());
//							}
//						}
						
//					}
					
	//			}else{
					/***********
					 * Over kids
					
					int px1 = 0;
					int px2 = 0;
					int x1 = 0;
					int x2 = 0;
					for (Integer kid : linkageMap.keySet()){
						px1 = Math.max(0, kid+dcp[0]);
						px2 = Math.min(currchromBinSize, kid+dcp[DCPsize-1]);
						x1 = (int) (1+(radius[n]-(kid-px1)));
						x2 = (int) (DCPsize-1 -(radius[n]-px2-kid));
						System.out.println("lenght of px is "+(px2-px1)+" length of x is "+(x2-x1));
						double maxScore = 0;
						int maxIndex = 0;
						*/
	//				}
				}
								
				//for each scaleNum, add the parents to the segmentationTree
				currScale.put(n, GvParents.keySet());
	
			}//end of scale space iteration
			
			for (Integer scale : currScale.keySet()){
				System.out.println("current scale is: "+scale);
				Set<Integer> nodesSet = currScale.get(scale);
				System.out.println("current nodeset size is: "+nodesSet.size());
				for (Integer node : nodesSet)
					System.out.println(node);
			}

			
	//		segmentationTree.put(currChrom, currScale);

			
			
	//		GaussianBlur = null;
			
		}// end of chromosome iteration		
		manager.close();
	}

		
	public static void main(String[] args) {
		
		/***
		 * Need to specify --tag3ext & --binwidth
		 ***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		SEEDConfig sconf = new SEEDConfig(gconf, args);
		DifferentialMSR profile = new DifferentialMSR (gconf, econf, sconf);	
		profile.buildMSR();
	}
	
}
