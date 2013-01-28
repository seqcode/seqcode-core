package edu.psu.compbio.seqcode.projects.multigps.mixturemodel;

import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

public class EMState{
	protected ExperimentManager manager;
	protected Config config;
	public double[][][] rBind; //Binding component responsibilities
	public double[][] rNoise;  //Noise component responsibilities
    public double[][] piBind;  // pi : emission probabilities for binding components
    public double[] piNoise;   // pi : emission probabilities for noise components
    public int[][] mu;         // mu : positions of the binding components

    public double LL;
    public int numNZComponent, numComponent;
    
    public EMState(ExperimentManager man, Config con, int totalComp, int[] hitNums){
    	manager=man;
    	config=con;
    	int numC = man.getNumConditions();
    	numComponent = totalComp;
    	mu = new int[numC][totalComp];
    	piBind = new double[numC][totalComp];
    	piNoise = new double[numC];
    	rBind = new double[numC][totalComp][];
    	rNoise = new double[numC][];
    	for(int c=0; c<numC; c++){
    		rNoise[c] = new double[hitNums[c]];
    		for(int x=0; x<totalComp; x++){
    			rBind[c][x] = new double[hitNums[c]];
    		}
    	}
    }
    // BIC=LL-#param/2*ln(n)
    // # param: Each component has 2 parameters, mixing prob and position, thus "*2";
    // "-1" comes from the fact that total mix prob sum to 1.
    // for multi-condition, # of beta variables is (mixture.numConditions-1)*numComponents
    // n: is the number of data point, i.e. the base positions.
    double BIC(double n){
        return LL - (numNZComponent*2-1 + (manager.getNumConditions()-1)*numNZComponent )/2*Math.log(n);
    }
    public String toString(){
        String str =  String.format("%.3f\t%d", LL, numNZComponent);
        for(int c=0; c<manager.getNumConditions(); c++){
	        str = str+"\n\t";
			for(int j=0; j<mu[c].length; j++){if(piBind[c][j]>0){
				str = str+piBind[c][j]+"\t";
			}}
		}
        return str;
    }
    public boolean equivalent(EMState s){
    	boolean numCompEqual = numNZComponent==s.numNZComponent;
    	boolean compPosEqual=true;
    	if(numCompEqual){
    		for(int c=0; c<manager.getNumConditions(); c++)
    			for(int j=0; j<mu[c].length; j++){if(piBind[c][j]>0){
    				compPosEqual = compPosEqual && (mu[c][j] == s.mu[c][j]);
    			}}
    	}else{
    		compPosEqual=false;
    	}
    	boolean piBindEquivalent=true;
    	for(int c=0; c<manager.getNumConditions(); c++)
			for(int j=0; j<piBind[c].length; j++){if(piBind[c][j]>0){
				piBindEquivalent = piBindEquivalent && (Math.abs(piBind[c][j]-s.piBind[c][j])<config.EM_STATE_EQUIV_THRES);
			}}
    	boolean rBindEquivalent=true;
    	for(int c=0; c<manager.getNumConditions(); c++)
    		for(int j=0; j<piBind[c].length; j++){if(piBind[c][j]>0){
    			for(int x=0; x<rBind[c][j].length; x++){
    				rBindEquivalent = rBindEquivalent && (Math.abs(rBind[c][j][x]-s.rBind[c][j][x])<config.EM_STATE_EQUIV_THRES);
    			}
			}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
    
    public void copy(EMState s){
    	int numC = manager.getNumConditions();
    	for(int c=0; c<numC; c++){
    		piNoise[c] = s.piNoise[c];
    		for(int x=0; x<s.rNoise[c].length; x++)
    			rNoise[c][x] = s.rNoise[c][x];
    		
    		for(int j=0; j<numComponent; j++){
    			piBind[c][j] = s.piBind[c][j];
    			mu[c][j] = s.mu[c][j];
    			for(int x=0; x<s.rBind[c][j].length; x++){
    				rBind[c][j][x] = s.rBind[c][j][x];
    			}
    		}
    	}
    }
}
