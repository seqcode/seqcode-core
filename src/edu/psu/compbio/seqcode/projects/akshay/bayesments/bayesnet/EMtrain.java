package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import java.util.Random;

import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentFeature;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;

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
	protected double[][][] Qijk;
	protected int numChromStates;
	protected int numFacBindingStates;
	protected int N; // number of training examples
	protected int C; // number of chromatin conditions
	protected int F; // number of factor conditions (alomost always 1)
	protected boolean finishedTraining;  // to know if the current state is trained or not
	protected boolean plot;
	
	public EMtrain(Config config, GenomicLocations trainingData, ExperimentManager manager) {
		this.config = config;
		this.trainingData = trainingData;
		this.plot = config.doEMplot();
		
		//Initializing the model
		initializeEM(manager);
		finishedTraining=false;
	}
	
	
	private void initializeEM(ExperimentManager manager){
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
		
		//Initialization from emperical means
		
		for(int c=0; c< C; c++){
			float[] observedValues = new float[N];
			for(int i=0; i<N; i++){
				observedValues[i] = Xc[i][c];
			}
			double[] means = this.getEmpMeanValues(numChromStates, observedValues);
			for(int j=0; j<numChromStates; j++){
				MUc[j][c] = means[j];
			}
		}
		
		for(int f=0; f< F; f++){
			float[] observedValues = new float[N];
			for(int i=0; i<N; i++){
				observedValues[i] = Xf[i][f];
			}
			double[] means = this.getEmpMeanValues(numFacBindingStates, observedValues);
			for(int k=0; k<numFacBindingStates; k++){
				MUf[k][f] = means[k];
			}
		}
		
		// random initialization
		//for(int i=0; i<numChromStates; i++){
		//	MUc[i] = this.getRandomList(C, false);
		//}
		//for(int i=0; i<numFacBindingStates; i++){
		//	MUf[i] = this.getRandomList(F, false);
		//}
		
		//Printing Mu's for debugging
		BayesmentsSandbox.printArray(MUc, "MUc", "MUc", manager);
		BayesmentsSandbox.printArray(MUf, "MUf", "MUf", manager);
		
	
		//Initializing sigma's
		SIGMAc = new double[numChromStates][C];
		SIGMAf = new double[numFacBindingStates][F];
		
		// Initializing for emperical data
		for(int c=0; c< C; c++){
			double[] observedValues = new double[N];
			for(int i=0; i<N; i++){
				observedValues[i] = Xc[i][c];
			}
			double min = observedValues[this.getMinindex(observedValues)];
			double max = observedValues[this.getMaxindex(observedValues)];
			System.out.println("Max "+ Double.toString(max));
			System.out.println("Min "+ Double.toString(min));
			for(int j=0; j<numChromStates; j++){
				SIGMAc[j][c] = max-min;
			}
		}
		
		for(int f=0; f< F; f++){
			double[] observedValues = new double[N];
			for(int i=0; i<N; i++){
				observedValues[i] = Xf[i][f];
			}
			double min = observedValues[this.getMinindex(observedValues)];
			double max = observedValues[this.getMaxindex(observedValues)];
			for(int k=0; k<numFacBindingStates; k++){
				SIGMAf[k][f] = max-min;
			}
		}
		
		
		//random initialization
		//for(int i=0; i<numChromStates; i++){
		//	SIGMAc[i] = this.getRandomList(C, false);
		//}
		//for(int i=0; i<numFacBindingStates; i++){
		//	SIGMAf[i] = this.getRandomList(F, false);
		//}
		
		//printing SIGMA's for debugging
		BayesmentsSandbox.printArray(SIGMAc, "SIGMAc", "SIGMAc", manager);
		BayesmentsSandbox.printArray(SIGMAf, "SIGMAf", "SIGMAf", manager);
		
		// Initializing Bjk
		Bjk = new double[numChromStates][numFacBindingStates];
		for(int i=0; i<numChromStates; i++){
			Bjk[i] = this.getRandomList(numFacBindingStates, true);
		}
		
		//printing Bjk for debugging
		BayesmentsSandbox.printArray(Bjk, "chrom_state", "factor_State", manager);
		
		// Initializing PIj
		PIj = new double[numChromStates];
		//PIj = this.getRandomList(numChromStates, true);
		PIj = this.getUniformList(numChromStates);
		//printing PIj for debugging
		BayesmentsSandbox.printArray(PIj, "chrom_state");
		
		//Initializing Qijk
		Qijk = new double[N][numChromStates][numFacBindingStates];
		
	}
	
	
	//Generates an array of random positive doubles that sum up to 1 
	private double[] getRandomList(int n, boolean prob){
		double[] ret = new double[n];
		double sum=0.0;
		Random ran = new Random();
		for(int i=0; i<n; i++){
			ret[i] = ran.nextDouble()+(double) ran.nextInt(90);
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
	
	private double[] getEmpMeanValues(int n, float[] observedValues){
		double[] ret =  new double[n];
		double mean = 0.0;
		double std=0.0;
		for(int i=0; i< observedValues.length; i++){
			mean = mean + observedValues[i];
		}
		mean = mean/(double) observedValues.length;
		
		for(int i=0; i<observedValues.length; i++){
			std = std + Math.pow(observedValues[i]-mean, 2.0);
		}
		std = std/(double) observedValues.length;
		
		
		std = Math.sqrt(std);
		Random rn =new Random();
		int range = (int) (0.4*std);
		
		range = range ==0 ? 1 : range;
		for(int i=0; i<n ; i++){
			double random  = rn.nextInt(range)+mean-0.2*std;
			ret[i] = random;
		}
		return ret;
	}
	
	private double[] getUniformList(int n){
		double[] ret =  new double[n];
		double value = 1/n;
		for(int i=0; i<n; i++){
			ret[i] = 1/(double)n;
		}
		return ret;
	}
	
	private int getMinindex(double[] list){
		double val = 100000.0;
		int ret=0;
		for(int i=0; i< list.length; i++){
			if(val > list[i]){
				ret = i;
				val = list[i];
			}
		}
		return ret;
	}
	
	private int getMaxindex(double[] list){
		double val = -100000.0;
		int ret=0;
		for(int i=0; i< list.length; i++){
			if(list[i] > val){
				ret =i;
				val = list[i];
			}
		}
		return ret;
	}
	
	public void runEM(){
		int itrs = config.getNumItrs();
		double[][][] trainMUc = new double[itrs+1][numChromStates][C]; //initial random params plus itrs 
		double[][][] trainMUf = new double[itrs+1][numFacBindingStates][F];
		double[][][] trainSIGMAc = new double[itrs+1][numChromStates][C];
		double[][][] trainSIGMAf = new double[itrs+1][numFacBindingStates][F];
		double[][] trainPIj = new double[itrs+1][numChromStates];
		double[][][] trainBjk = new double[itrs+1][numChromStates][numFacBindingStates];
		
		for(int t=0; t<itrs; t++){
			System.out.println("Iteration no: "+ Integer.toString(t));
			if(t==0){      //copy the initial set of random params
				trainMUc[0] = MUc;
				trainMUf[0] = MUf;
				trainSIGMAc[0] = SIGMAc;
				trainSIGMAf[0] = SIGMAf;
				trainPIj[0] = PIj;
				trainBjk[0] = Bjk;
			}
			
			executeEStep();  //E-Step
			executeMStep(); //M-Step
			
			// copy the updated params
			for(int j=0; j<numChromStates; j++){
				for(int c=0; c<C; c++){
					trainMUc[t+1][j][c]= MUc[j][c];
				}
			}
			
			for(int k=0; k<numFacBindingStates; k++){
				for(int f=0; f<F; f++){
					trainMUf[t+1][k][f]= MUf[k][f];
				}
			}
			
			for(int j=0; j<numChromStates; j++){
				for(int c=0; c<C; c++){
					trainSIGMAc[t+1][j][c]= SIGMAc[j][c];
				}
			}
			
			for(int k=0; k<numFacBindingStates; k++){
				for(int f=0; f<F; f++){
					trainSIGMAf[t+1][k][f]= SIGMAf[k][f];
				}
			}
			
			for(int j=0; j<numChromStates; j++){
				trainPIj[t+1][j] = PIj[j];
			}
			
			for(int j=0; j< numChromStates; j++){
				for(int k=0; k<numFacBindingStates; k++){
					trainBjk[t+1][j][k] = Bjk[j][k];
				}
			}
			
		}
		
		this.finishedTraining=true;
		
		// plot if asked for
		if(plot){
			EMplotter ep = new EMplotter(config, trainMUc, trainMUf, trainSIGMAc, trainSIGMAf, trainPIj, trainBjk, C, F);
		}
		
		
		
	}
	
	private void executeEStep(){
		//E-Step
		double den[]= new double[N];
		for(int i=0; i<N; i++){    // over the training examples
			for(int j=0; j<numChromStates; j++){   // over the chromatin states
				for(int k=0; k< numFacBindingStates; k++){ //over factor binding states
					double chromGausssianProd=0.0;
					double facGaussianProd = 0.0;
					for(int c=0; c<C; c++){
						NormalDistribution gaussian = new NormalDistribution(MUc[j][c],Math.pow(SIGMAc[j][c], 2.0));
						chromGausssianProd = (c==0 ? gaussian.calcProbability((double) Xc[i][c]): chromGausssianProd* gaussian.calcProbability((double) Xc[i][c]));
					}
					for(int f=0; f< F; f++){
						NormalDistribution gaussian = new NormalDistribution(MUf[k][f],Math.pow(SIGMAf[k][f], 2.0));
						facGaussianProd = (f == 0 ? gaussian.calcProbability((double) Xf[i][f]): facGaussianProd* gaussian.calcProbability((double) Xf[i][f]));
					}
					chromGausssianProd = ( Double.isNaN(chromGausssianProd)) ? 0.0 : chromGausssianProd;
					facGaussianProd = (Double.isNaN(facGaussianProd)) ? 0.0: facGaussianProd;
					Qijk[i][j][k] = PIj[j]*chromGausssianProd*Bjk[j][k]*facGaussianProd;
					//debug lines
					//System.out.println("pie value: "+Double.toString(PIj[j]));
					//System.out.println("chromatin guassian product value: "+Double.toString(chromGausssianProd));
					//System.out.println("factor guassian product value: "+Double.toString(facGaussianProd));
					//System.out.println("Bjk value: "+Double.toString(Bjk[j][k]));
					//System.out.println("product value: "+Double.toString(Qijk[i][j][k]));
					den[i] = den[i]+Qijk[i][j][k];
				}
				//debug line
				//System.out.println("den for i: "+Double.toString(den[i]));
			}
			
		}
		
		for(int i=0; i<N; i++){
			for(int j=0; j<numChromStates; j++){
				for(int k=0; k<numFacBindingStates; k++){
					Qijk[i][j][k] = Qijk[i][j][k]/den[i];
					
					Qijk[i][j][k] = ( Double.isNaN(Qijk[i][j][k])) ? 0.0 : Qijk[i][j][k];
					//debug line
					//System.out.println("qijk values: "+Double.toString(Qijk[i][j][k]));
				}
			}
		}
		
	}
	
	private void executeMStep(){
		
		//-------------------------PIj update-----------------------------
		
		//compute
		double denPIj = 0.0;
		for(int j=0; j<numChromStates; j++){
			for(int i=0; i<N; i++){
				for(int k=0; k<numFacBindingStates; k++){
					PIj[j] = k==0 && i==0 ? Qijk[i][j][k] : PIj[j]+Qijk[i][j][k];
					denPIj = denPIj + Qijk[i][j][k];
				}
			}
			
			//debug lines
			//System.out.println("PIj values: "+Double.toString(PIj[j]));
		}
		
		//normalize
		for(int j=0; j<numChromStates; j++){
			PIj[j] = PIj[j]/denPIj;
		}
		
		//making sure PI-j for ant state does not go to zero
		
		for(int j=0; j<numChromStates; j++){
			if(PIj[j] == 0){
				PIj[j] = 0.001;
			}
		}
			//re-normalize
		denPIj = 0.0;
		for(int j=0; j< numChromStates; j++){
			denPIj = denPIj + PIj[j];
		}
		for(int j=0; j<numChromStates; j++){
			PIj[j] = PIj[j]/denPIj;
		}
		
		//debug line
		//System.out.println("denom for PIj");
		//System.out.println(Double.toString(denPIj));
		
		//debug lines
		//System.out.println("PI-j");
		//BayesmentsSandbox.printArray(PIj, "chrom");
		
		//-----------------------MUc update------------------------------------
		
		//compute
		double[][] denMUc=new double[numChromStates][C];
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				for(int i=0; i<N; i++){
					for(int k=0; k<numFacBindingStates; k++){
						MUc[j][c] = k==0 && i==0 ? Qijk[i][j][k]*Xc[i][c] : MUc[j][c]+ Qijk[i][j][k]*Xc[i][c];
						denMUc[j][c] = denMUc[j][c]+Qijk[i][j][k];
					
					}
				}
				
				//debug line
				//System.out.println("muec num values: "+Double.toString(MUc[j][c]));
				//System.out.println("muec din vakues: "+ Double.toString(denMUc[j][c]));
			}
		}
		
		//dividing by denominator 
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				
				//debug line
				if(denMUc[j][c] == 0){System.out.println("Chromatin state "+Integer.toString(j)+" Condition "+Integer.toString(c));}
				MUc[j][c] = MUc[j][c]/denMUc[j][c];
			}
		}
		
		
		//-----------------------MUf update --------------------------------------
		
		//compute
		double[][] denMUf = new double[numFacBindingStates][F];
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				for(int i=0; i<N; i++){
					for(int j=0; j<numChromStates; j++){
						MUf[k][f] = j==0 && i==0? Qijk[i][j][k]*Xf[i][f]: MUf[k][f]+Qijk[i][j][k]*Xf[i][f];
						denMUf[k][f] = denMUf[k][f]+Qijk[i][j][k];
					}
				}
				//debug line
				//System.out.println("muef num values: "+Double.toString(MUf[k][f]));
				//System.out.println("muef din vakues: "+ Double.toString(denMUf[k][f]));
			}
		}
		
		//dividing by denominator
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				MUf[k][f] = MUf[k][f]/denMUf[k][f];
			}
		}
		
		//debug lines
		//BayesmentsSandbox.printArray(MUc, "chrom_state", "Condition");
		//BayesmentsSandbox.printArray(MUf, "factor_state", "Condition");
		
		//--------------------SIGMAc update --------------------------------------
		
		//compute
		double[][] denSIGMAc = new double[numChromStates][C];
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				for(int i=0; i<N; i++){
					for(int k=0; k<numFacBindingStates; k++){
						SIGMAc[j][c] = k==0 && i==0? Qijk[i][j][k]*Math.pow(((double)Xc[i][c] - MUc[j][c]) , 2.0): SIGMAc[j][c]+Qijk[i][j][k]*Math.pow(((double)Xc[i][c] - MUc[j][c]) , 2.0);
						denSIGMAc[j][c] = denSIGMAc[j][c]+Qijk[i][j][k];
					}
				}
				//debug line
				//System.out.println("sigmac num values: "+Double.toString(SIGMAc[j][c]));
				//System.out.println("sigmac din vakues: "+ Double.toString(denSIGMAc[j][c]));
			}
		}
		
		//dividing by denominator and taking square root
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				
				
				//debug line
				if(denSIGMAc[j][c] == 0){System.out.println("Chromatin state "+Integer.toString(j)+" Condition "+Integer.toString(c));}
				SIGMAc[j][c] = Math.sqrt(SIGMAc[j][c]/denSIGMAc[j][c]);
			}
		}
		
		//------------------SIGMAf update------------------------------------------
		
		//compute
		double[][] denSIGMAf = new double[numFacBindingStates][F];
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				for(int i=0; i<N; i++){
					for(int j=0; j< numChromStates; j++){
						SIGMAf[k][f] = j==0 && i==0? Qijk[i][j][k]*Math.pow(((double)Xf[i][f] - MUf[k][f]) , 2.0) : SIGMAf[k][f]+Qijk[i][j][k]*Math.pow(((double)Xf[i][f] - MUf[k][f]) , 2.0);
						denSIGMAf[k][f] = denSIGMAf[k][f]+Qijk[i][j][k];
					}
				}
				
				//debug line
				//System.out.println("sigmaf num values: "+Double.toString(SIGMAf[k][f]));
				//System.out.println("sigmaf din vakues: "+ Double.toString(denMUc[k][f]));
			}
		}
		
		//dividing by denominator and taking square root
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				SIGMAf[k][f] = Math.sqrt(SIGMAf[k][f]/denSIGMAf[k][f]);
			}
		}
		
		//printing SIGMA's for debugging
		//BayesmentsSandbox.printArray(SIGMAc, "chrom_state", "Condition");
		//BayesmentsSandbox.printArray(SIGMAf, "factor_state", "Condition");
		
		//---------------------Bjk update -------------------------------------------
		
		//compute
		double[] denBjk = new double[numFacBindingStates];
		for(int k=0; k<numFacBindingStates; k++){
			for(int j=0; j< numChromStates; j++){
				for(int i=0; i<N; i++){
					Bjk[j][k] = i==0 ? Qijk[i][j][k] : Bjk[j][k]+Qijk[i][j][k];
					denBjk[k] = denBjk[k]+Qijk[i][j][k];
				}
			}
		}
		
		//normalize
		for(int k=0; k<numFacBindingStates; k++){
			for(int j=0; j< numChromStates; j++){
				Bjk[j][k] = Bjk[j][k]/denBjk[k];
			}
		}
		//printing Bjk for debugging
		//BayesmentsSandbox.printArray(Bjk, "chrom_state", "factor_State");
	}
	
	//Accessors
	
	public boolean isTrainined(){return this.finishedTraining;}
	public double[] getPIj(){return this.PIj;}
	public double[][] getMUc(){return this.MUc;}
	public double[][] getMUf(){return this.MUf;}
	public double[][] getSIGMAc(){return this.SIGMAc;}
	public double[][] getSIGMAf(){return this.SIGMAf;}
	public double[][] getBjk(){return this.Bjk;}
	public GenomicLocations getTrainingData(){return this.trainingData;}
	
	// main method is only for testing puposers
	
	public static void main(String[] args){
		//double[] test = getRandomList(3,true);
		//System.out.println(test[0]);
		//System.out.println(test[1]);
		//System.out.println(test[2]);
	}
	
}
