package edu.psu.compbio.seqcode.projects.naomi;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.distribution.NormalDistribution;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;

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

	public DifferentialMSR(GenomeConfig gcon, ExptConfig econ){	
		gconfig = gcon;
		econfig = econ;
	}
	 
	public void buildMSR(){
		
		/*********************
		 * Gaussian scale space and window parameters	
		 */
		// arbitrary number of scale
		int numScale= 10;
		double DELTA_TAU = 0.5*Math.log(2);	
		// I have to determine P_MIN value carefully because P_MIN will substantially affect Gaussian window size
		double P_MIN = Math.pow(10,-3);
		double K_MIN = 1/Math.sqrt(1-Math.exp(-2*DELTA_TAU));	
		double K_N = Math.ceil(K_MIN);
		
		/*********************
		 * Starting nodes
		 */
		HashSet <Integer> kids = new HashSet<Integer>();
		
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
		
		Iterator<Region> chroms = new ChromosomeGenerator<Genome>().execute(genome);
		//iterating each chromosome (each chromosome is a region).
		while (chroms.hasNext()) {
			
			Region currChrom = chroms.next();			
			int currchromSize= currChrom.getWidth()+1;
			
			//primitive matrix to store signal and the subsequent convolved signals
			//its index correspond to the coordinates
			float[][] GaussianBlur = new float[currchromSize][2];
			for (int i = 0;i<currchromSize;i++){
				for (int j = 0; j<2;j++)
					GaussianBlur[i][j] = 0;
			}
			
			//get StrandedBaseCount list for each chromosome
			Map<Sample, List<StrandedBaseCount>> sampleCountsMap = new HashMap<Sample, List<StrandedBaseCount>>();
			for (Sample sample : manager.getSamples())
				sampleCountsMap.put(sample,sample.getBases(currChrom)); 
			
			//StrandedBasedCount object contains positive and negative strand separately
			//store all base counts indexed by positions at column[1]
			//Also detect kids (starting nodes)
			
			for (Sample sample : manager.getSamples()){		
				List<StrandedBaseCount> currentCounts = sampleCountsMap.get(sample);
				for (StrandedBaseCount hits: currentCounts)
					GaussianBlur[hits.getCoordinate()][1]+=hits.getCount();
				currentCounts = null;
			}
			
			//adding starting nodes; to qualify for the starting nodes the signal intensity needs to be different from the subsequent signal intensity
			//also add the starting and end positions in the kids at start and end positions  
			kids.add(0);
			for (int i = 1 ; i< GaussianBlur.length-1; i++){
				if (GaussianBlur[i] != GaussianBlur[i+1])
						kids.add(i);
			}
			kids.add(GaussianBlur.length-1);

			for (int n = 2;n < numScale+1; n++){
				
				/*********************
				 * Gaussian scale space 	
				 */
				double polyCoeffi[] = new double [currchromSize];
				//first copy from column[1] to column[0];this procedure need to be repeated for each iteration of scale
				//also copy from column[1] to array to store polynomial coefficient
				for (int i = 0 ; i<currchromSize; i++){
					GaussianBlur[i][0]=GaussianBlur[i][1];
					polyCoeffi[i]=GaussianBlur[i][1];
				}
				//sigma calculation
				sigma[n] = Math.exp((n-1)*DELTA_TAU);
				// create normal distribution with mean zero and sigma[n]
				NormalDistribution normDistribution = new NormalDistribution(0.00,sigma[n]);
				//take inverse CDF based on the normal distribution using probability
				double inverseCDF = normDistribution.inverseCumulativeProbability(P_MIN);
				int windowSize = (int) (-Math.round(inverseCDF)*2+1);
				double X[]=new double[windowSize];
				for (int i = 0; i<windowSize;i++)
					X[i] = Math.round(inverseCDF)+i;
			
				//window calculation based on Gaussian(normal) distribution with sigma, mean=zero,x=X[i]			
				double window[] = new double[windowSize];
				double windowSum = 0;
				for (int i = 0;i<windowSize;i++){
					window[i] = normDistribution.density(X[i]);
					windowSum=windowSum+window[i];
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
				
				//copy Gaussian blur results to the column[1]
				// I should check to make sure that it's not off by 1
				for (int i = 0; i<currchromSize;i++){
					if (currchromSize % 2 ==0 && coefficients.length/2 == 1)
						GaussianBlur[i][1]=(float) coefficients[polyMid-currchromSize/2+i+1];
					else
						GaussianBlur[i][1]=(float) coefficients[polyMid-currchromSize/2+i];
				}	
				
				/***************
				 * Search Volume	
				 */
				double tempRadius;
				if (n==2){
					tempRadius = sigma[n];
				}else{
					tempRadius = Math.sqrt(Math.pow(sigma[n],2)-Math.pow(sigma[n-1], 2));
				}
				radius[n] = Math.ceil(K_MIN*tempRadius);
					
				int DCPsize = (int) (Math.round(radius[n])*2+1);
				double dcp, distance;
				double distanceFactor[] = new double[DCPsize];
				double denom = -2*(Math.pow(sigma[n], 2)-Math.pow(sigma[n-1],2));
				for (int i = 0; i<DCPsize;i++){
					dcp = -radius[n];
					// applying equation 7 in Vincken(1997)
					distance=Math.exp(dcp/denom)/Math.exp(0.5*Math.pow(sigma[n],2)/denom);
					
					//applying equation 8 in Vincken (1997) 
					if (Math.abs(dcp) > 0.5*sigma[n]){
						distanceFactor[i]= distance;
					}else{
						distanceFactor[i] = 1.0000;
					}
					dcp = dcp+1;
				}
				
				/***************
				 * Linkage Loop	
				 */
				int kids_len = kids.size();
				double parents [] = new double[kids_len];
				for (int i = 0; i<kids_len;i++)
					parents[i] = 0;
				
				if (DCPsize < 50){				
					/***********
					 * Over window
					 */
					
				}
				
				
 
					
					
				
							
				
				
			}
		}

	}
		
	public static void main(String[] args) {
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
	
	}
	
}
