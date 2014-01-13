package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import java.util.Random;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentFeature;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class EMtrain {

	protected Config config;
	protected GenomicLocations trainingData;
	protected float[][] Xc;
	protected float[][] Xf;
	protected double[][] MUc;
	protected double[][] MUf;
	protected double[][] SIGMAc;
	protected double[][] SIGMAf;
	protected double[][] Bjk;
	protected double[] PIj;
	protected int numChromStates;
	protected int numFacBindingStates;
	protected int N; // number of training examples
	protected int C; // number of chromatin conditions
	protected int F; // number of factor conditions (alomost always 1)
	protected boolean plot;
	
	public EMtrain(Config config, GenomicLocations trainingData) {
		this.config = config;
		this.trainingData = trainingData;
		this.plot = config.doEMplot();
		
		//Initializing the model
		initializeEM();
	}
	
	
	public void initializeEM(){
		// getting N, M, C, P and F
		N=this.trainingData.getNumTrainingExamples();
		C=this.trainingData.getNumChromatinCons();
		F=this.trainingData.getNumFacCons();
		numChromStates = config.getNumChrmStates();
		numFacBindingStates = config.getNumFacStates();
		
		//Initializing and loading X's
		this.Xc = new float[N][C];
		this.Xf = new float[N][F];
		this.Xc = this.trainingData.getChromatinCounts();
		this.Xf = this.trainingData.getFactorCounts();
		
		//Initializing mu's
		MUc = new double[numChromStates][C];
		MUf = new double[numFacBindingStates][F];
		for(int i=0; i<numChromStates; i++){
			MUc[i] = this.getRandomList(C, false);
		}
		for(int i=0; i<numFacBindingStates; i++){
			MUf[i] = this.getRandomList(F, false);
		}
		
		//Initializing sigma's
		SIGMAc = new double[numChromStates][C];
		SIGMAf = new double[numFacBindingStates][F];
		for(int i=0; i<numChromStates; i++){
			SIGMAc[i] = this.getRandomList(C, false);
		}
		for(int i=0; i<numFacBindingStates; i++){
			SIGMAf[i] = this.getRandomList(C, false);
		}
		
		// Initializing Bjk
		Bjk = new double[numChromStates][numFacBindingStates];
		for(int i=0; i<numChromStates; i++){
			Bjk[i] = this.getRandomList(numFacBindingStates, true);
		}
		
		// Initializing PIj
		PIj = new double[numChromStates];
		PIj = this.getRandomList(numChromStates, true);
		
	}
	
	
	//Generates an array of random positive doubles that sum up to 1 
	private double[] getRandomList(int n, boolean prob){
		double[] ret = new double[n];
		double sum=0.0;
		Random ran = new Random();
		for(int i=0; i<n; i++){
			ret[i] = ran.nextDouble()+(double) ran.nextInt();
			sum = sum + ret[i];
		}
		if(prob){
			for(int i=0; i<n; i++){
				ret[i] = ret[i]/sum;
			}
			return ret;
		}else{
			return ret;
		}
	}
	
	
}
