package edu.psu.compbio.seqcode.projects.multigps.mixturemodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.BackgroundCollection;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;
import edu.psu.compbio.seqcode.projects.multigps.utilities.EMStepPlotter;

/**
 * BindingMLAssignment: Maximum likelihood assignment of reads to a configuration of binding components.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class BindingMLAssignment {

	protected ExperimentManager manager;
	protected Config config;
	protected List<BindingComponent> components;
	protected List<NoiseComponent> noise;
	protected int numComponents;  //Assumes the same number of active+inactive components in each condition
	protected int numConditions;
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds; //Background models per condition
	//	EM VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[][]   sigHitCounts;	// Hit weights
	protected int[][]      sigHitPos;		// Hit positions
	protected boolean[][]  sigHitPlusStr;	// Hit positive strand boolean
	protected int[]		   sigHitNum;		// Number of hits in each condition 
	protected int[][]      sigRepIndices;  // Index of replicate for the hit
	protected double[][]   ctrlHitCounts;	// Hit weights
	protected int[][]      ctrlHitPos;		// Hit positions
	protected boolean[][]  ctrlHitPlusStr;	// Hit positive strand boolean
	protected int[]		   ctrlHitNum;		// Number of hits in each condition 
	protected int[][]      ctrlRepIndices;  // Index of replicate for the hit
	protected double[][][] h; 			// H function (binding component probability per read)
	protected double[][]   n; 			// N function (noise component probability per read)
	protected double[][][] rBindSig;		// Binding component responsibilities (signal reads)
	protected double[][]   rNoiseSig;		// Noise component responsibilities (signal reads)
	protected double[][][] rBindCtrl;		// Binding component responsibilities (control reads)
	protected double[][]   rNoiseCtrl;		// Noise component responsibilities (control reads)
	protected double[][]   pi;			// pi : emission probabilities for binding components
	protected double[]     piNoise;		// pi : emission probabilities for noise components (fixed)
	protected int[][]      mu;			// mu : positions of the binding components
	protected double[][]   compLL;		//Log-likelihood for each component in each condition
	protected double[][][] lastRBind;	//Last responsibilities (monitor convergence)
	protected double[][]   lastPi;		//Last Pi (monitor convergence)
	protected int[][]      lastMu;		//Last positions (monitor convergence)
	protected double[][]   tmp_pi;			// pi used in ML calc
	protected double[][][] tmp_h;			// h used in ML calc
	protected double[][][] tmp_rBindSig;	// rBindSig used in ML calc
	protected double[][]   tmp_rNoiseSig;	// rNoiseSig used in ML calc
	protected double[]	   tmp_piNoise; 	// piNoise used in ML calc
	protected BindingModel[] bindingModels; //Array of binding models for convenience
	protected double[]	   sigRepHitCountTotals; //Hit count totals counted by replicate (for convenience)
	protected double[]	uniformRepHitCountTotals; //Hit count totals by replicate if signal read counts were distributed uniformly (used only if there is no control) 
	protected double numPotentialRegions;
	
	/**
	 * Constructor
	 * @param c
	 * @param eMan
	 */
	public BindingMLAssignment(Config c, ExperimentManager eMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg){
		config=c;
		manager = eMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;
	}
	
	/**
     * ML assignment
     *
     * Takes as input a SINGLE list of binding components (i.e. a single configuraton)
     * Returns lists of binding events indexed by condition. 
     * Pi values are calculated from the signal hits and then those same components are directly applied to the control hits. 
     *
     * Almost purely matrix/array operations.
     */
    public List<BindingEvent>  assign(List<List<StrandedBaseCount>> signals,
    								  List<List<StrandedBaseCount>> controls,
    								  Region w, 
    								  List<NoiseComponent> noise,
    								  List<BindingComponent> comps, 
    								  int numComp){
    	components = comps;
        this.noise = noise;
        numComponents = numComp;
        //Matrix initializations
        sigHitCounts= new double[numConditions][];	// Hit weights
    	sigHitPos= new int[numConditions][];			// Hit positions
    	sigHitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	sigHitNum = new int[numConditions];			// Number of hits in each condition
    	sigRepIndices= new int[numConditions][];	    // Index of replicate for the hit
    	ctrlHitCounts= new double[numConditions][];	// Hit weights
    	ctrlHitPos= new int[numConditions][];			// Hit positions
    	ctrlHitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	ctrlHitNum = new int[numConditions];			// Number of hits in each condition
    	ctrlRepIndices= new int[numConditions][];	    // Index of replicate for the hit
    	h= new double[numConditions][][]; 			// H function (binding component probability per read)
    	n= new double[numConditions][]; 			// N function (noise component probability per read)
    	rBindSig= new double[numConditions][][];		// Binding component responsibilities (signal reads)
    	rNoiseSig= new double[numConditions][];		// Noise component responsibilities (signal reads)
    	rBindCtrl= new double[numConditions][][];		// Binding component responsibilities (control reads)
    	rNoiseCtrl= new double[numConditions][];		// Noise component responsibilities (control reads)
    	pi = new double[numConditions][numComponents];	// pi : emission probabilities for binding components
    	piNoise = new double[numConditions];		// pi : emission probabilities for noise components (fixed)
    	mu = new int[numConditions][numComponents];// mu : positions of the binding components
    	compLL = new double [numConditions][numComponents];		//Log-likelihood for each component in each condition
        bindingModels = new BindingModel[manager.getExperimentSet().getReplicates().size()]; //Array of bindingModels for convenience
        sigRepHitCountTotals = new double[manager.getExperimentSet().getReplicates().size()]; //Hit count totals counted by replicate (for convenience)
        uniformRepHitCountTotals = new double[manager.getExperimentSet().getReplicates().size()]; //Hit count totals by replicate if reads were distributed uniformly 
        //Monitor state convergence using the following last variables
        lastRBind = new double[numConditions][][];
        lastPi = new double[numConditions][numComponents];
        lastMu = new int[numConditions][numComponents];
        //Temporary variables
        tmp_pi = new double[numConditions][numComponents];			// pi used in ML calc
    	tmp_rBindSig= new double[numConditions][][];	// rBindSig used in ML calc
    	tmp_rNoiseSig= new double[numConditions][];	// rNoiseSig used in ML calc
    	tmp_piNoise = new double[numConditions]; 	// piNoise used in ML calc
    	tmp_h= new double[numConditions][][]; 			// H function used in ML calc
        
        
        //Initializing data structures
        for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
        	int c = cond.getIndex();
        	
        	//Add bindingModels to array
        	for(ControlledExperiment rep : cond.getReplicates())
        		bindingModels[rep.getIndex()] = rep.getBindingModel();
        	
        	//Load Reads (merge from all replicates)
        	List<StrandedBaseCount> sigBases = new ArrayList<StrandedBaseCount>();
        	List<StrandedBaseCount> ctrlBases = new ArrayList<StrandedBaseCount>();
        	for(ControlledExperiment rep : cond.getReplicates()){
        		sigBases.addAll(signals.get(rep.getIndex()));
        		if(controls.get(rep.getIndex())!=null)
        			ctrlBases.addAll(controls.get(rep.getIndex()));
        	}
        	sigHitNum[c] = sigBases.size();
        	ctrlHitNum[c] = ctrlBases.size();
        	
        	//Count total weights for convenience
        	for(ControlledExperiment rep : cond.getReplicates()){
        		sigRepHitCountTotals[rep.getIndex()]=0;
        		for(StrandedBaseCount s : signals.get(rep.getIndex()))
        			sigRepHitCountTotals[rep.getIndex()]+=s.getCount();
        		uniformRepHitCountTotals[rep.getIndex()] = ((rep.getNoiseCount()/config.getMappableGenomeLength())*(double)w.getWidth())/rep.getControlScaling(); //Normalizing by control scaling is a hack - usually control scaling will be 1 when the replicate has no control... however, it is not 1 for SES. 
        	}
        	
        	//Load replicate index for each read
        	sigRepIndices[c] = new int[sigHitNum[c]];
        	ctrlRepIndices[c] = ctrlHitNum[c]==0 ? null : new int[ctrlHitNum[c]];
        	int ys=0, yc=0, z=0;
        	for(ControlledExperiment rep : cond.getReplicates()){
        		z=0;
        		while(z<signals.get(rep.getIndex()).size()){
        			sigRepIndices[c][ys] = rep.getIndex();
        			z++; ys++;
        		}
        		z=0;
        		while(z<controls.get(rep.getIndex()).size()){
        			ctrlRepIndices[c][yc] = rep.getIndex();
        			z++; yc++;
        		}
        	}
            
            //Load signal read info
            sigHitCounts[c]= new double[sigHitNum[c]];
            sigHitPos[c]= new int[sigHitNum[c]];
            sigHitPlusStr[c]= new boolean[sigHitNum[c]];
            for(int i=0;i<sigHitNum[c];i++){
            	sigHitPos[c][i] = sigBases.get(i).getCoordinate();
            	sigHitPlusStr[c][i] = sigBases.get(i).getStrand() == '+';
            	sigHitCounts[c][i]=sigBases.get(i).getCount();
            }
            
            //Load control read info
            if(ctrlHitNum[c]>0){
	            ctrlHitCounts[c]= new double[ctrlHitNum[c]];
	            ctrlHitPos[c]= new int[ctrlHitNum[c]];
	            ctrlHitPlusStr[c]= new boolean[ctrlHitNum[c]];
	            for(int i=0;i<ctrlHitNum[c];i++){
	            	ctrlHitPos[c][i] = ctrlBases.get(i).getCoordinate();
	            	ctrlHitPlusStr[c][i] = ctrlBases.get(i).getStrand() == '+';
	            	ctrlHitCounts[c][i]=ctrlBases.get(i).getCount();
	            }
            }

            //Load pi for binding components
            for(int j=0;j<numComp;j++){
                BindingComponent comp = components.get(j);
                pi[c][j]= comp.getPi(); 
            }
            //Load pi for noise components
            piNoise[c]=noise.get(c).getPi();
            
            //Load binding component positions
            for(int j=0;j<numComp;j++)
            	mu[c][j] = components.get(j).getPosition();
    		
            //Initialize responsibility functions
            double[][] hc= new double[numComp][sigHitNum[c]];
            double[][] thc= new double[numComp][sigHitNum[c]];
            double[] nc = new double[sigHitNum[c]];
            for(int i=0;i<sigHitNum[c];i++){
            	for(int j=0;j<numComp;j++){
            		int dist = sigHitPlusStr[c][i] ? sigHitPos[c][i]-mu[c][j]: mu[c][j]-sigHitPos[c][i];
                    hc[j][i] = bindingModels[sigRepIndices[c][i]].probability(dist);
                    thc[j][i] = bindingModels[sigRepIndices[c][i]].probability(dist);
                }
            	nc[i] = noise.get(c).scorePosition(sigHitPos[c][i],sigRepIndices[c][i]);
            }
            h[c] = hc;
            n[c] = nc;
            tmp_h[c] = thc;
            
            rBindSig[c]  = new double[numComp][sigHitNum[c]];
    		rNoiseSig[c] = new double[sigHitNum[c]];
    		lastRBind[c] = new double[numComp][sigHitNum[c]];
    		tmp_rBindSig[c]  = new double[numComp][sigHitNum[c]];
    		tmp_rNoiseSig[c] = new double[sigHitNum[c]];
        }//End of data structure initialization
        
        
        //////////
        // Run ML steps
        //////////
        ML(w);
        
        
        //////////
        // Assign ML result to BindingEvents
        //////////
        List<BindingEvent> events = new ArrayList<BindingEvent>();
	   	for(int j=0;j<numComp;j++){ 
	   		BindingEvent event = new BindingEvent(components.get(j).getCoord(), w);
	   		
    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
    			ArrayList<Sample> controlsSeen = new ArrayList<Sample>();
    			boolean uniformBackAdded=false;
    			double condSigResp = 0.0, condCtrlResp = 0.0;
    			int c = cond.getIndex();
    			for(ControlledExperiment rep : cond.getReplicates()){
    				int r = rep.getIndex();
	    			double repSigResp = 0.0, repCtrlResp = 0.0; 
	    			if(pi[c][j]>0){
	    				double scount=0;
			            for(int i=0;i<sigHitNum[c];i++)
			            	if(sigRepIndices[c][i]==r)
			            		scount += sigHitCounts[c][i]*rBindSig[c][j][i];
			            repSigResp+=scount;
			            condSigResp+=scount;
			            
			            double ccount=0;
			            if(rep.hasControl()){
			            	for(int i=0;i<ctrlHitNum[c];i++)
			            		if(ctrlRepIndices[c][i]==r)
			            			ccount += ctrlHitCounts[c][i]*rBindCtrl[c][j][i];
			            	repCtrlResp+=ccount;
				            if(!controlsSeen.contains(rep.getControl()))
				            	condCtrlResp+=ccount;
				            controlsSeen.add(rep.getControl());
			            }else{  //If there is no control channel, assign pseudo-control counts as if the noise reads in the IP channel were distributed perfectly uniformly
			            	repCtrlResp = uniformRepHitCountTotals[rep.getIndex()]*pi[c][j];
			            	if(!uniformBackAdded)
			            		condCtrlResp+=repCtrlResp;
			            	uniformBackAdded=true;
			            }
	    			}
		            event.setRepSigHits(rep, repSigResp);
		            event.setRepCtrlHits(rep, repCtrlResp);
    			}
    			event.setCondSigHits(cond, condSigResp);
	            event.setCondCtrlHits(cond, condCtrlResp);
	            if(config.CALC_COMP_LL)
	            	event.setLLd(cond, compLL[c][j]);
	    	}
    		if(nonZeroInAny(event))
    			events.add(event);
        }        
        return events;
    }//end of EMTrain method

    private boolean nonZeroInAny(BindingEvent ev){
    	boolean nonZero=false;
    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
    		for(ControlledExperiment rep : cond.getReplicates()){
    			nonZero = nonZero || ev.getRepSigHits(rep)>=1;
    		}
    	}
    	return nonZero;
    }

    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void ML (Region currRegion) {
        int numComp = numComponents;
        double [][] totalRespSig = new double[numConditions][];
        double [][] totalRespCtrl = new double[numConditions][];
        
        //Initialize responsibilities
        for(int c=0; c<numConditions; c++){
    		totalRespSig[c] = new double[sigHitNum[c]];
    		for(int i=0;i<sigHitNum[c];i++)
                totalRespSig[c][i] = 0;
    		if(ctrlHitNum[c]>0){
	    		rBindCtrl[c] = new double[numComp][ctrlHitNum[c]];
	    		rNoiseCtrl[c]= new double[ctrlHitNum[c]];
	    		totalRespCtrl[c] = new double[ctrlHitNum[c]];
	            for(int i=0;i<ctrlHitNum[c];i++)
	                totalRespCtrl[c][i] = 0;
    		}
    	}        
    	
    	////////////////////////////
        //Run ML -- this should only need one or two rounds
    	////////////////////////////
        for(int t=0; t<config.EM_ML_ITER ; t++){ //System.out.println(t); 
        	
    		////////
    		//E-step
    		////////
    		for(int c=0; c<numConditions; c++){ int numBases = sigHitNum[c];
        		//Recompute h function, given binding component positions (n function is constant because noise model doesn't move)
        		for(int i=0;i<numBases;i++)
                	for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
                    	int dist = sigHitPlusStr[c][i] ? sigHitPos[c][i]-mu[c][j]: mu[c][j]-sigHitPos[c][i];
                    	h[c][j][i] = bindingModels[sigRepIndices[c][i]].probability(dist);
                	}}
        		//Compute responsibilities
    			for(int i=0;i<numBases;i++)
                    totalRespSig[c][i] = 0;
        		for(int i=0;i<numBases;i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				rBindSig[c][j][i] = h[c][j][i]*pi[c][j];
        				totalRespSig[c][i] +=rBindSig[c][j][i]; 
        			}}
        			rNoiseSig[c][i] = n[c][i] * piNoise[c];
        			totalRespSig[c][i] +=rNoiseSig[c][i];
        		}
        		//Normalize responsibilities
        		for(int i=0;i<numBases;i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				rBindSig[c][j][i]/=totalRespSig[c][i];
        			}}
        			rNoiseSig[c][i]/=totalRespSig[c][i];
        		}
    		}
    		        		
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
    		for(int c=0; c<numConditions; c++){ int numBases = sigHitNum[c];	
        		//Maximize pi
        		double[] sumR=new double[numComponents];
        		for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        			for(int i=0;i<numBases;i++)
        				sumR[j] += rBindSig[c][j][i]*sigHitCounts[c][i];
                }}
                
        		// No components to be eliminated in ML, update pi(j)
        		for(int j=0;j<numComp;j++){ 
        			pi[c][j]=Math.max(0, sumR[j]); 
        		}
                
                //Normalize pi (accounting for piNoise)
                double totalPi=0;
                for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
            		totalPi+=pi[c][j];
            	}}
                for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
            		if(totalPi>0)
            			pi[c][j]=pi[c][j]/(totalPi/(1-piNoise[c]));
            	}}
        	}
        	
        	//Non-zero components count
        	int nonZeroComps=0;
        	for(int c=0; c<numConditions; c++)
        		for(int j=0;j<numComp;j++)
        			if(pi[c][j]>0.0)
        				nonZeroComps++;
        	
        	////////////
        	//Check Stopping condition
        	////////////	
            if (nonZeroComps>0 && (t==0 || !lastEquivToCurr())){
            	copyStateToLast();
                continue;
            }else{
            	copyStateToLast();
            	break;
            }
        } //LOOP: Run ML while not converged
        //Base log-likelihood calculation
        double[] baseLL =new double[numConditions];
        for(int c=0; c<numConditions; c++){
        	baseLL[c]=0;
    		int numBases = sigHitNum[c];
    		for(int i=0;i<numBases;i++){
    			// for each read, each event will give a conditional prob or bg prob
                double j_sum=0;
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0.0){
    				j_sum += Math.log(rBindSig[c][j][i])/config.LOG2;
                }}
    			j_sum += Math.log(rNoiseSig[c][i])/config.LOG2;
                baseLL[c] += j_sum*sigHitCounts[c][i];                        
            }
        }
        
        //ML assignment of signal reads to components is finished
        //Assign control reads with converged pi values here
        for(int c=0; c<numConditions; c++){ int numBases = ctrlHitNum[c];
        	double[][] hCtrl= new double[numComp][numBases];
            double[] nCtrl = new double[numBases];

			//Recompute h & n functions for control reads, given binding component positions 
			for(int i=0;i<numBases;i++){
	        	for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
	            	int dist = ctrlHitPlusStr[c][i] ? ctrlHitPos[c][i]-mu[c][j]: mu[c][j]-ctrlHitPos[c][i];
	            	hCtrl[j][i] = bindingModels[ctrlRepIndices[c][i]].probability(dist);
	        	}}
	        	nCtrl[i] = noise.get(c).scorePosition(ctrlHitPos[c][i], ctrlRepIndices[c][i]);
			}
			//Compute responsibilities
			for(int i=0;i<numBases;i++)
	            totalRespCtrl[c][i] = 0;
			for(int i=0;i<numBases;i++){
				for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
					rBindCtrl[c][j][i] = hCtrl[j][i]*pi[c][j];
					totalRespCtrl[c][i] +=rBindCtrl[c][j][i]; 
				}}
				rNoiseCtrl[c][i] = nCtrl[i] * piNoise[c];
				totalRespCtrl[c][i] +=rNoiseCtrl[c][i];
			}
			//Normalize responsibilities
			for(int i=0;i<numBases;i++){
				for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
					rBindCtrl[c][j][i]/=totalRespCtrl[c][i];
				}}
				rNoiseCtrl[c][i]/=totalRespCtrl[c][i];
			}
		}
        
        if(config.CALC_COMP_LL){
	        /////////////////////////////////////////////////////////////////////////////////
	        //Finally, run ML on alternate models, where each component is eliminated in turn
	        /////////////////////////////////////////////////////////////////////////////////
	        for(int elim=0;elim<numComp;elim++){ 
	        	//If this component is active in any condition
	        	boolean active=false;
	        	for(int c=0; c<numConditions; c++){if(pi[c][elim]>0){ active=true;}}
	        	if(active){
	        		for(int t=0; t<=2 ; t++){  
	        			//Initialize temporary variables
	        			for(int c=0; c<numConditions; c++){
	        				for(int j=0;j<numComp;j++)
	        					if(j==elim)
	        						tmp_pi[c][j]=0;
	        					else
	        						tmp_pi[c][j]= pi[c][j];
	        				if(piNoise[c]+pi[c][elim]==1.0)//Single component case
	        					tmp_piNoise[c]=1.0;
	        				else{
	        					tmp_piNoise[c] = piNoise[c];
	        					double totalPi=0;
	                            for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                        		totalPi+=tmp_pi[c][j];
	                        	}}
	                            for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                        		if(totalPi>0)
	                        			tmp_pi[c][j]=tmp_pi[c][j]/(totalPi/(1-tmp_piNoise[c]));
	                        	}}	
	        				}
	        			}
	        			
	                	////////
	            		//E-step
	            		////////
	            		for(int c=0; c<numConditions; c++){ int numBases = sigHitNum[c];
	                		//Recompute h function, given binding component positions (n function is constant because noise model doesn't move)
	                		for(int i=0;i<numBases;i++)
	                        	for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                            	int dist = sigHitPlusStr[c][i] ? sigHitPos[c][i]-mu[c][j]: mu[c][j]-sigHitPos[c][i];
	                            	tmp_h[c][j][i] = bindingModels[sigRepIndices[c][i]].probability(dist);
	                        	}}
	                		//Compute responsibilities
	            			for(int i=0;i<numBases;i++)
	                            totalRespSig[c][i] = 0;
	                		for(int i=0;i<numBases;i++){
	                			for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                				tmp_rBindSig[c][j][i] = tmp_h[c][j][i]*tmp_pi[c][j];
	                				totalRespSig[c][i] +=tmp_rBindSig[c][j][i]; 
	                			}}
	                			tmp_rNoiseSig[c][i] = n[c][i] * tmp_piNoise[c];
	                			totalRespSig[c][i] +=rNoiseSig[c][i];
	                		}
	                		//Normalize responsibilities
	                		for(int i=0;i<numBases;i++){
	                			for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                				tmp_rBindSig[c][j][i]/=totalRespSig[c][i];
	                			}}
	                			tmp_rNoiseSig[c][i]/=totalRespSig[c][i];
	                		}
	            		}
	            		        
	            		/////////////////////
	            		//M-step: maximize pi
	            		/////////////////////
	            		for(int c=0; c<numConditions; c++){ int numBases = sigHitNum[c];
	            			if(tmp_piNoise[c]<1.0){//No need for maximization in single component cases
	                		//Maximize pi
	                		double[] sumR=new double[numComponents];
	                		for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
	                			for(int i=0;i<numBases;i++)
	                				sumR[j] += tmp_rBindSig[c][j][i]*sigHitCounts[c][i];
	                        }}
	                        
	                		// No components to be eliminated in ML, update pi(j)
	                		for(int j=0;j<numComp;j++){ 
	                			tmp_pi[c][j]=Math.max(0, sumR[j]); 
	                		}
	                        
	                        //Normalize pi (accounting for piNoise)
	                        double totalPi=0;
	                        for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                    		totalPi+=tmp_pi[c][j];
	                    	}}
	                        for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0){
	                    		if(totalPi>0)
	                    			tmp_pi[c][j]=tmp_pi[c][j]/(totalPi/(1-tmp_piNoise[c]));
	                    	}}
	            			}
	                	}
	        		}
	        		// log-likelihood calculation
	                for(int c=0; c<numConditions; c++){
	                	compLL[c][elim] =-2*baseLL[c];
	            		int numBases = sigHitNum[c];
	            		for(int i=0;i<numBases;i++){
	            			// for each read, each event will give a conditional prob or bg prob
	                        double j_sum=0;
	            			for(int j=0;j<numComp;j++){ if(tmp_pi[c][j]>0.0){
	            				j_sum += Math.log(tmp_rBindSig[c][j][i])/config.LOG2;
	                        }}
	            			j_sum += Math.log(tmp_rNoiseSig[c][i])/config.LOG2;
	                        compLL[c][elim] += 2*j_sum*sigHitCounts[c][i];                        
	                    }
	                }
	        	}
	        }
        }        
    }//end of EM_MAP method
 
    /**
     * Copy current variables to last variables (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     */
    private void copyStateToLast(){
    	int numC = manager.getNumConditions();
    	for(int c=0; c<numC; c++){
    		for(int j=0; j<numComponents; j++){
    			lastPi[c][j] = pi[c][j];
    			lastMu[c][j] = mu[c][j];
    			for(int x=0; x<rBindSig[c][j].length; x++){
    				lastRBind[c][j][x] = rBindSig[c][j][x];
    			}
    		}
    	}
    }
    
    /**
     * Compare last variables to current (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     * @return
     */
    private boolean lastEquivToCurr(){
    	int numC = manager.getNumConditions();
    	int currNZ=0, lastNZ=0;
    	for(int c=0; c<numConditions; c++)
    		for(int j=0;j<numComponents;j++){
    			if(pi[c][j]>0)
    				currNZ++;
    			if(lastPi[c][j]>0)
    				lastNZ++;
    		}
    	boolean numCompEqual = currNZ==lastNZ;
    	boolean compPosEqual=true;
    	if(numCompEqual){
    		for(int c=0; c<numC; c++)
    			for(int j=0; j<mu[c].length; j++){if(pi[c][j]>0){
    				compPosEqual = compPosEqual && (mu[c][j] == lastMu[c][j]);
    			}}
    	}else{
    		compPosEqual=false;
    	}
    	boolean piBindEquivalent=true;
    	for(int c=0; c<numC; c++)
			for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
				piBindEquivalent = piBindEquivalent && (Math.abs(pi[c][j]-lastPi[c][j])<config.EM_STATE_EQUIV_THRES);
			}}
    	boolean rBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
    			for(int x=0; x<rBindSig[c][j].length; x++){
    				rBindEquivalent = rBindEquivalent && (Math.abs(rBindSig[c][j][x]-lastRBind[c][j][x])<config.EM_STATE_EQUIV_THRES);
    			}
			}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
}
