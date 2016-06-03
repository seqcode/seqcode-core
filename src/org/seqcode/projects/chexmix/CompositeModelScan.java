package org.seqcode.projects.chexmix;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;


/**
 * CompositeModelScan: scan a set of regions for the best XL matches
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CompositeModelScan {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ChExMixConfig config;
	protected ExperimentManager manager;
	protected ProteinDNAInteractionModel model; //The model to scan with
	protected WindowedTagDistributions tagDistributions; //The windows to scan
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[][]   hitCounts;	// Hit weights
	protected int[][]      hitPos;		// Hit positions
	protected boolean[][]  hitPlusStr;	// Hit positive strand boolean
	protected int[]		   hitNum;		// Number of hits in each condition 
	protected double[][][] hAll;		// H function values for all positions in the current window (precomputed)
	protected double[][][] h; 			// H function (binding component probability per read)
	protected double[][][] r;		// Binding component responsibilities
	protected double[]   pi;			// pi : emission probabilities for binding components
	protected int[]      mu;			// mu : positions of the binding components
	protected int numConditions;
	protected int numComponents;
	protected int scanWindow;
	
	public CompositeModelScan(ProteinDNAInteractionModel model, int scanWindow, GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon, ExperimentManager eMan){
		gconfig = gcon;
		econfig = econ;
		config = ccon;
		manager = eMan;
		this.model = model;
		numConditions = manager.getNumConditions();
		numComponents = model.getNumComponents();
		this.scanWindow = scanWindow;
	}
	
	//Accessors
	public ProteinDNAInteractionModel getModel(){return model;}
	
	
	/**
	 * Run the scanner  
	 * 
	 */
	public List<StrandedPoint> scan(WindowedTagDistributions tagDists){		
		List<StrandedPoint> scanMaxPoints = new ArrayList<StrandedPoint>();
		this.tagDistributions = tagDists;
	
		//Matrix initializations
        hitCounts= new double[numConditions][];	// Hit weights
    	hitPos= new int[numConditions][];			// Hit positions
    	hitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	hitNum = new int[numConditions];			// Number of hits in each condition
    	h= new double[numConditions][][]; 			// H function (binding component probability per read)
    	r= new double[numConditions][][];		// Binding component responsibilities
    	pi = new double[numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numComponents]; // mu : positions of the binding components (fixed across conditions)
        
        //Initializing data structures that will be used in every scan
        for(int j=0;j<numComponents;j++){
        	//Load pi for binding components
        	CompositeModelComponent comp = model.getAllComponents().get(j);
            pi[j]= comp.getPi();
        }
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
        	//Number of unique positions in the composite distribution (plus & minus)
        	hitNum[c]=tagDistributions.getWinSize()*2;
            //Count matrices
            int[] posc= new int[hitNum[c]];
            boolean[] plusc= new boolean[hitNum[c]];
            for(int i=0;i<tagDistributions.getWinSize();i++){ //Watson
            	posc[i] = i;
            	plusc[i] = true;
            }
            for(int i=0;i<tagDistributions.getWinSize();i++){ //Crick
            	posc[i+tagDistributions.getWinSize()] = i;
            	plusc[i+tagDistributions.getWinSize()] = false;                    
            }
            hitPos[c] = posc;
            hitPlusStr[c] = plusc;
            r[c] = new double[numComponents][hitNum[c]];
        }
           
        //Iterate through each region
        for(Region reg : tagDistributions.getRegions()){
        	double bestScore = -Double.MAX_VALUE;
        	int bestPos = -1;
        	int bestStrand=0;

        	//Scan each strand
        	for(int strand=0; strand<2; strand++){
	        	//Initialize this region's data
	        	for(ExperimentCondition cond : manager.getConditions()){
	            	int c = cond.getIndex();
	            	//Count matrices
	                double[] countc= new double[hitNum[c]];
	                if(strand==1){
		                for(int i=0;i<tagDistributions.getWinSize();i++)
		                    countc[i]=tagDistributions.getRegionPlus(reg, cond)[i];
		                for(int i=0;i<tagDistributions.getWinSize();i++)
		                    countc[i+tagDistributions.getWinSize()]=tagDistributions.getRegionMinus(reg, cond)[i];
	                }else{
	                	for(int i=0;i<tagDistributions.getWinSize();i++)
		                    countc[i]=tagDistributions.getRegionMinus(reg, cond)[tagDistributions.getWinSize()-i-1];
		                for(int i=0;i<tagDistributions.getWinSize();i++)
		                    countc[i+tagDistributions.getWinSize()]=tagDistributions.getRegionPlus(reg, cond)[tagDistributions.getWinSize()-i-1];
	                }
	                hitCounts[c]=countc;
	        	}
	        	
	        	//Scan through window
	        	for(int scanOffset=0; scanOffset<scanWindow; scanOffset++){
	        		//Load binding component positions
	        		for(int j=0;j<numComponents;j++)
	        			mu[j] = model.getAllComponents().get(j).getPosition()+scanOffset; //needs to be re-initialized for every window shift
	        		
	        		
	        		//Initialize responsibility functions
	                for(ExperimentCondition cond : manager.getConditions()){
	                	int c = cond.getIndex();
	                	double[][] hc= new double[numComponents][hitNum[c]];
	                    for(int i=0;i<hitNum[c];i++){
	                    	for(int j=0;j<numComponents;j++){ if(pi[j]>0){
	                    		if(hitPos[c][i]>=scanOffset && hitPos[c][i]<=(tagDistributions.getWinSize()-scanWindow+scanOffset)){
	                    			int dist = hitPos[c][i]-mu[j];
	                    			hc[j][i] = model.getAllComponents().get(j).getTagDistribution().probability(dist, hitPlusStr[c][i]);
	                    		}else{
	                    			hc[j][i] = 0;
	                    		}
	                        }}
	                    }
	                    h[c] = hc;        	
	                }
	                
					//////////////////////////////////////////////////////////////////////////////////////////
					//Run E step once
					//////////////////////////////////////////////////////////////////////////////////////////
					double sumResp = 0;
	                for(int c=0; c<numConditions; c++){ 
						//Compute responsibilities
						//for(int i=0;i<hitNum[c];i++)
						//	totalResp[c][i] = 0;
						for(int i=0;i<hitNum[c];i++){
							//for(int j=0;j<numComponents;j++){ if(pi[j]>0){
							for(CompositeModelComponent xlComp : model.getXLComponents()){
								int j=xlComp.getIndex(); 
								if(pi[j]>0){
									r[c][j][i] = h[c][j][i]*pi[j];
									sumResp+=r[c][j][i]*hitCounts[c][i];
									//totalResp[c][i] +=sumResp;
								}
							}
						}
						//Normalize responsibilities
						//for(int i=0;i<hitNum[c];i++){
						//	for(int j=0;j<numComponents;j++){ if(pi[j]>0){
						//		r[c][j][i]/=totalResp[c][i];
						//	}}
						//}
					}
					
					//Test for best match!
					if(sumResp>bestScore){
						bestScore = sumResp;
						bestPos = scanOffset;
						bestStrand = strand;
					}
	        	}
        	}
        	StrandedPoint bestPoint = bestStrand==0 ? 
        			new StrandedPoint(reg.getGenome(), reg.getChrom(), (reg.getStart()+bestPos+model.getCenterOffset()), '+') :
        			new StrandedPoint(reg.getGenome(), reg.getChrom(), (reg.getEnd()-bestPos-model.getCenterOffset()), '-');
        	scanMaxPoints.add(bestPoint);
        }
		return scanMaxPoints;
	}
	

}
