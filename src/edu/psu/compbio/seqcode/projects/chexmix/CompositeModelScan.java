package edu.psu.compbio.seqcode.projects.chexmix;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;

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
	
	public CompositeModelScan(ProteinDNAInteractionModel model, GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon, ExperimentManager eMan){
		gconfig = gcon;
		econfig = econ;
		config = ccon;
		manager = eMan;
		this.model = model;
		numConditions = manager.getNumConditions();
		numComponents = model.getNumComponents();
	}
	
	//Accessors
	public ProteinDNAInteractionModel getModel(){return model;}
	
	
	/**
	 * Run the scanner  
	 * 
	 */
	public void scan(WindowedTagDistributions tagDists){		
		this.tagDistributions = tagDists;
	/*	
		//Matrix initializations
        hitCounts= new double[numConditions][];	// Hit weights
    	hitPos= new int[numConditions][];			// Hit positions
    	hitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	hitNum = new int[numConditions];			// Number of hits in each condition
    	h= new double[numConditions][][]; 			// H function (binding component probability per read)
    	r= new double[numConditions][][];		// Binding component responsibilities
    	pi = new double[numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numComponents]; // mu : positions of the binding components (fixed across conditions)
        
        //Initializing data structures
        for(int j=0;j<numComponents;j++){
        	//Load pi for binding components
        	CompositeModelComponent comp = model.getAllComponents().get(j);
            pi[j]= comp.getPi();
            //Load binding component positions
        	mu[j] = model.getAllComponents().get(j).getPosition()+scanOffset; //needs to be re-initialized for every window shift
        }
        
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
        	
        	//Number of unique positions in the composite distribution (plus & minus)
        	hitNum[c]=tagDistributions.getWinSize()*2;

            //Load read info
            double[] countc= new double[hitNum[c]];
            int[] posc= new int[hitNum[c]];
            boolean[] plusc= new boolean[hitNum[c]];
            for(int i=0;i<tagDistributions.getWinSize();i++){ //Watson
            	posc[i] = i;
            	plusc[i] = true;
                countc[i]=tagDistributions.getCompositeWatson(cond)[i]; //change to site-specific
            }
            for(int i=0;i<tagDistributions.getWinSize();i++){ //Crick
            	posc[i+tagDistributions.getWinSize()] = i;
            	plusc[i+tagDistributions.getWinSize()] = false;
                countc[i+tagDistributions.getWinSize()]=tagDistributions.getCompositeCrick(cond)[i];
            }
            hitPos[c] = posc;
            hitCounts[c]=countc;
            hitPlusStr[c] = plusc;
        	
            //Initialize responsibility functions
            double[][] hc= new double[numComponents][hitNum[c]];
            for(int i=0;i<hitNum[c];i++){
            	for(int j=0;j<numComponents;j++){
            		int dist = hitPos[c][i]-mu[j];
                    hc[j][i] = model.getAllComponents().get(j).getTagDistribution().probability(dist, hitPlusStr[c][i]);
                }
            }
            h[c] = hc;
            
            r[c] = new double[numComponents][hitNum[c]];
    		        	
        }
        
        double [][] totalResp = new double[numConditions][];
        
        //Initialize responsibilities
        for(int c=0; c<numConditions; c++){
    		int numBases = hitNum[c];
    		totalResp[c] = new double[numBases];
            for(int i=0;i<numBases;i++)
                totalResp[c][i] = 0;
    	}        
    	
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run E step once
    	//////////////////////////////////////////////////////////////////////////////////////////
        for(int c=0; c<numConditions; c++){ 
    		//Recompute h function, given binding component positions
    		for(int i=0;i<hitNum[c];i++)
            	for(int j=0;j<numComponents;j++){ if(pi[j]>0){
            		int dist = hitPos[c][i]-mu[j];
            		h[c][j][i] = model.getAllComponents().get(j).getTagDistribution().probability(dist, hitPlusStr[c][i]);
            	}}
    		//Compute responsibilities
			for(int i=0;i<hitNum[c];i++)
                totalResp[c][i] = 0;
    		for(int i=0;i<hitNum[c];i++){
    			for(int j=0;j<numComponents;j++){ if(pi[j]>0){
    				r[c][j][i] = h[c][j][i]*pi[j];
    				totalResp[c][i] +=r[c][j][i]; 
    			}}
    		}
    		//Normalize responsibilities
    		for(int i=0;i<hitNum[c];i++){
    			for(int j=0;j<numComponents;j++){ if(pi[j]>0){
    				r[c][j][i]/=totalResp[c][i];
    			}}
    		}
		}
		*/
	}
	

}
