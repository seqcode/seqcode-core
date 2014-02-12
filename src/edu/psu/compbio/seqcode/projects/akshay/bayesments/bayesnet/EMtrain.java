package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import java.util.Random;

import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.Sequences;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;

/**
 * EMtrain class the initializes and trains the parameters of the Bayesian network using the EM framework
 * @author akshaykakumanu
 *
 */

public class EMtrain {

	protected Config config;
	protected GenomicLocations trainingData;
	
	protected Sequences seqs;
	protected double[][] Xs;
	protected double[][] MUs;
	protected double[][] SIGMAs;
	
	protected boolean seqState;
	
	// 2-d array of chromatin counts, with rows as locations and columns as conditions
	protected float[][] Xc;
	//2-d array of factor countts, with rows as locations and columns as conditions
	protected float[][] Xf;
	//2-d array of means, with rows as chromatin states and columns as conditions
	protected double[][] MUc;
	//2-d array of means, with rows as factor stats and columns as conditions
	protected double[][] MUf;
	//2-d array of variance, with rows as chromatin states and columns as conditions
	protected double[][] SIGMAc;
	//2-d array of variance, with rows as factor stats and columns as conditions
	protected double[][] SIGMAf;
	//2-d array of transition probabilities, with rows as chromatin states and colums as factor states
	protected double[][] Bjk;
	//1-d array of probabilities for different chromatin states
	protected double[] PIj;
	//Expectation of unobserved variables for the E step in EM algorithm
	protected double[][][] Qijk;
	//no of chromatin states
	protected int numChromStates;
	//number of factor states
	protected int numFacBindingStates;
	protected int N; // number of training examples
	protected int C; // number of chromatin conditions
	protected int F; // number of factor conditions (alomost always 1)
	protected int M; // number of motifs
	// flag to turn on plotting of the parameters as the function of iterations step in EM
	protected boolean plot;
	
	double[][][] trainMUs;
	double[][][] trainSIGMAs;
	
	protected int itr_no=0;
	protected int total_itrs=0;
	
	/**
	 * When initializing this class object, it is always in seq off mode
	 * @param config
	 * @param trainingData
	 * @param manager
	 */
	public EMtrain(Config config, GenomicLocations trainingData, ExperimentManager manager, int total_itrs) {
		this.config = config;
		this.trainingData = trainingData;
		this.plot = config.doEMplot();
		
		//Initializing the model
		initializeEM(manager);
		this.total_itrs = total_itrs;
	}
	
	/**
	 * Method that initializes all the parameters for the Bayesian network
	 * @param manager
	 */
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
		
		//Printing the initial Mu's
		BayesmentsSandbox.printArray(MUc, "MUc", "MUc", manager);
		BayesmentsSandbox.printArray(MUf, "MUf", "MUf", manager);
		
		//Initializing sigma's
		SIGMAc = new double[numChromStates][C];
		SIGMAf = new double[numFacBindingStates][F];
		
		// Initializing from emperical data... (Max-Min) 
		for(int c=0; c< C; c++){
			double[] observedValues = new double[N];
			for(int i=0; i<N; i++){
				observedValues[i] = Xc[i][c];
			}
			double min = observedValues[this.getMinindex(observedValues)];
			double max = observedValues[this.getMaxindex(observedValues)];
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
		
		//Printing the initial SIGMA's
		BayesmentsSandbox.printArray(SIGMAc, "SIGMAc", "SIGMAc", manager);
		BayesmentsSandbox.printArray(SIGMAf, "SIGMAf", "SIGMAf", manager);
		
		// Initializing Bjk .. Using random initialization 
		Bjk = new double[numChromStates][numFacBindingStates];
		for(int i=0; i<numChromStates; i++){
			Bjk[i] = this.getRandomList(numFacBindingStates, true);
		}
		
		//printing the initial Bjk's
		BayesmentsSandbox.printArray(Bjk, "chrom_state", "factor_State", manager);
		
		// Initializing PIj ... Using Uniform initialization
		PIj = new double[numChromStates];
		PIj = this.getUniformList(numChromStates);
		
		//Printing the initial PIj's
		BayesmentsSandbox.printArray(PIj, "chrom_state");
		
		//Initializing the dimensions of Qijk's
		Qijk = new double[N][numChromStates][numFacBindingStates];
		
	}
	
	private void initializeSeqParams(){
		
		this.MUs = new double[numChromStates][M];
		this.SIGMAs = new double[numChromStates][M];
		
		for(int m=0; m< M; m++){
			float[] observedScores = new float[N];
			for(int i=0; i<N; i++){
				observedScores[i] = (float) Xs[i][m];
			}
			double[] means = this.getEmpMeanValues(numChromStates,  observedScores);
			for(int j=0; j<numChromStates; j++){
				MUc[j][m] = means[j];
			}
		}
		
		for(int m=0; m< M; m++){
			double[] observedValues = new double[N];
			for(int i=0; i<N; i++){
				observedValues[i] = Xs[i][m];
			}
			double min = observedValues[this.getMinindex(observedValues)];
			double max = observedValues[this.getMaxindex(observedValues)];
			
			for(int j=0; j<numChromStates; j++){
				SIGMAc[j][m] = max-min;
			}
		}
	}
	
	
	/**
	 * Generates n random numbers that sum up to 1 if the boolean prob is TRUE 
	 * @param n
	 * @param prob
	 * @return
	 */
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
	
	/**
	 * Generates n values that fall in mu+-0.2*sigma of the the observedValues
	 * @param n
	 * @param observedValues
	 * @return
	 */
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
		double range = 0.4*std;
		
		for(int i=0; i<n ; i++){
			double random  = rn.nextDouble()*range+mean-0.2*std;
			ret[i] = random;
		}
		return ret;
	}
	
	/**
	 * Returns n equal numbers that add upto one
	 * @param n
	 * @return
	 */
	private double[] getUniformList(int n){
		double[] ret =  new double[n];
		double value = 1/n;
		for(int i=0; i<n; i++){
			ret[i] = 1/(double)n;
		}
		return ret;
	}
	
	/**
	 * Returns the index of the minimum number in the given list of doubles
	 * @param list
	 * @return
	 */
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
	
	/**
	 * returns the index of the maximum number in the given list of doubles
	 * @param list
	 * @return
	 */
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
	
	/**
	 * Runs the EM algorithm for a given number of iterations and also plots the parameters over the learning rounds if plot parameter is true
	 */
	public void runEM(int itrs){
		
		// Initializing the arrays to store the parameters for all training rounds
		double[][][] trainMUc = new double[this.total_itrs+1][numChromStates][C]; //initial random params plus itrs 
		double[][][] trainMUf = new double[this.total_itrs+1][numFacBindingStates][F];
		double[][][] trainSIGMAc = new double[this.total_itrs+1][numChromStates][C];
		double[][][] trainSIGMAf = new double[this.total_itrs+1][numFacBindingStates][F];
		double[][] trainPIj = new double[this.total_itrs+1][numChromStates];
		double[][][] trainBjk = new double[this.total_itrs+1][numChromStates][numFacBindingStates];
		
		
		for(int t=0; t<itrs; t++){ // training for the given number of iterations
			
			if(itr_no==0){      //Copy the initial set of random parameters
				trainMUc[0] = MUc;
				trainMUf[0] = MUf;
				trainSIGMAc[0] = SIGMAc;
				trainSIGMAf[0] = SIGMAf;
				trainPIj[0] = PIj;
				trainBjk[0] = Bjk;
				itr_no++;
			}
			if(this.seqState){
				if(t == 0){
					for(int p=0; p< this.itr_no-1 ; p++){
						for(int j=0; j<this.numChromStates; j++){
							for(int m=0; m<M; m++){
									trainMUs[p][j][m] =0.0; 
							}
						}
					}
					for(int j=0; j<this.numChromStates; j++){
						for(int m=0; m<M; m++){
							trainMUs[itr_no-1] = MUs;
							trainSIGMAs[itr_no-1] = SIGMAs;
						}
					}
				}
			}
			
			executeEStep();  //E-Step
			executeMStep(); //M-Step
			
			// Copy the updated parameters
			for(int j=0; j<numChromStates; j++){
				for(int c=0; c<C; c++){
					trainMUc[itr_no][j][c]= MUc[j][c];
				}
			}
			
			for(int k=0; k<numFacBindingStates; k++){
				for(int f=0; f<F; f++){
					trainMUf[itr_no][k][f]= MUf[k][f];
				}
			}
			
			for(int j=0; j<numChromStates; j++){
				for(int c=0; c<C; c++){
					trainSIGMAc[itr_no][j][c]= SIGMAc[j][c];
				}
			}
			
			for(int k=0; k<numFacBindingStates; k++){
				for(int f=0; f<F; f++){
					trainSIGMAf[itr_no][k][f]= SIGMAf[k][f];
				}
			}
			
			for(int j=0; j<numChromStates; j++){
				trainPIj[itr_no][j] = PIj[j];
			}
			
			for(int j=0; j< numChromStates; j++){
				for(int k=0; k<numFacBindingStates; k++){
					trainBjk[itr_no][j][k] = Bjk[j][k];
				}
			}
			
			if(this.seqState){
				for(int j=0; j<numChromStates; j++){
					for(int m=0; m<M; m++){
						trainMUs[itr_no][j][m]= MUs[j][m];
					}
				}
				for(int j=0; j<numChromStates; j++){
					for(int m=0; m<M; m++){
						trainSIGMAs[itr_no][j][m]= SIGMAs[j][m];
					}
				}
			}
			this.itr_no++;
		}
		
		// Plot if asked for
		if(plot){
			EMplotter ep = new EMplotter(config, trainMUc, trainMUf, trainSIGMAc, trainSIGMAf, trainPIj, trainBjk, C, F);
		}
		
		
		
	}
	
	/**
	 * Executes the E step. That is it calculates the Qijk's using the current set of parameters
	 */
	private void executeEStep(){
		
		double den[]= new double[N];
		
		//Calculate the numerator and the denominator for all the Qijk's
		for(int i=0; i<N; i++){    // over the training examples
			for(int j=0; j<numChromStates; j++){   // over the chromatin states
				for(int k=0; k< numFacBindingStates; k++){ //over factor binding states
					double chromGausssianProd=0.0;
					double facGaussianProd = 0.0;
					double seqGaussianProd = 0.0;
					for(int c=0; c<C; c++){
						NormalDistribution gaussian = new NormalDistribution(MUc[j][c],Math.pow(SIGMAc[j][c], 2.0));
						chromGausssianProd = (c==0 ? gaussian.calcProbability((double) Xc[i][c]): chromGausssianProd* gaussian.calcProbability((double) Xc[i][c]));
					}
					for(int f=0; f< F; f++){
						NormalDistribution gaussian = new NormalDistribution(MUf[k][f],Math.pow(SIGMAf[k][f], 2.0));
						facGaussianProd = (f == 0 ? gaussian.calcProbability((double) Xf[i][f]): facGaussianProd* gaussian.calcProbability((double) Xf[i][f]));
					}
					// Set the guassian products to 0 in case they are NaN
					chromGausssianProd = ( Double.isNaN(chromGausssianProd)) ? 0.0 : chromGausssianProd;
					facGaussianProd = (Double.isNaN(facGaussianProd)) ? 0.0: facGaussianProd;
					
					if(this.seqState){
						for(int m=0; m<M; m++){
							NormalDistribution gaussian = new NormalDistribution(MUs[j][m],Math.pow(SIGMAs[j][m], 2.0));
							seqGaussianProd = (m==0 ? gaussian.calcProbability((double) Xs[i][m]): seqGaussianProd* gaussian.calcProbability((double) Xs[i][m]));
						}
						seqGaussianProd = (Double.isNaN(seqGaussianProd)) ? 0.0: seqGaussianProd;
					}
					
					if(this.seqState){
						Qijk[i][j][k] = PIj[j]*chromGausssianProd*Bjk[j][k]*facGaussianProd*seqGaussianProd;
					}else{Qijk[i][j][k] = PIj[j]*chromGausssianProd*Bjk[j][k]*facGaussianProd;}
					
					
					den[i] = den[i]+Qijk[i][j][k];
				}
			}
			
		}
		
		//Normalize the numerator by dividing the Qijk's with the denominators
		for(int i=0; i<N; i++){
			for(int j=0; j<numChromStates; j++){
				for(int k=0; k<numFacBindingStates; k++){
					Qijk[i][j][k] = Qijk[i][j][k]/den[i];
					
					Qijk[i][j][k] = ( Double.isNaN(Qijk[i][j][k])) ? 0.0 : Qijk[i][j][k];
				}
			}
		}
		
	}
	
	/**
	 * Executes the M step. Updates all the parameters of the model
	 */
	private void executeMStep(){
		
		//-------------------------PIj update-----------------------------
		
		//Compute
		double denPIj = 0.0;
		for(int j=0; j<numChromStates; j++){
			for(int i=0; i<N; i++){
				for(int k=0; k<numFacBindingStates; k++){
					PIj[j] = k==0 && i==0 ? Qijk[i][j][k] : PIj[j]+Qijk[i][j][k];
					denPIj = denPIj + Qijk[i][j][k];
				}
			}
		}
		
		//Normalize
		for(int j=0; j<numChromStates; j++){
			PIj[j] = PIj[j]/denPIj;
		}
		
		//Making sure PI-j for any state does not go to zero
		for(int j=0; j<numChromStates; j++){
			if(PIj[j] == 0){
				PIj[j] = 0.001;
			}
		}
		
		//Re-normalize
		denPIj = 0.0;
		for(int j=0; j< numChromStates; j++){
			denPIj = denPIj + PIj[j];
		}
		for(int j=0; j<numChromStates; j++){
			PIj[j] = PIj[j]/denPIj;
		}
		
		//-----------------------MUc update------------------------------------
		
		//Compute
		double[][] denMUc=new double[numChromStates][C];
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				for(int i=0; i<N; i++){
					for(int k=0; k<numFacBindingStates; k++){
						MUc[j][c] = k==0 && i==0 ? Qijk[i][j][k]*Xc[i][c] : MUc[j][c]+ Qijk[i][j][k]*Xc[i][c];
						denMUc[j][c] = denMUc[j][c]+Qijk[i][j][k];
					
					}
				}
			}
		}
		
		//Normalize 
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				MUc[j][c] = MUc[j][c]/denMUc[j][c];
			}
		}
		
		
		//-----------------------MUf update --------------------------------------
		
		//Compute
		double[][] denMUf = new double[numFacBindingStates][F];
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				for(int i=0; i<N; i++){
					for(int j=0; j<numChromStates; j++){
						MUf[k][f] = j==0 && i==0? Qijk[i][j][k]*Xf[i][f]: MUf[k][f]+Qijk[i][j][k]*Xf[i][f];
						denMUf[k][f] = denMUf[k][f]+Qijk[i][j][k];
					}
				}
			}
		}
		
		//Normalize
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				MUf[k][f] = MUf[k][f]/denMUf[k][f];
			}
		}
		
		// -----------------------MUs update, if in seqState----------------------------
		
		//Compute
		if(this.seqState){
			double[][] denMUs = new double[numChromStates][M];
			for(int j=0; j<numChromStates; j++){
				for(int m=0; m<M; m++){
					for(int i=0; i<N; i++){
						for(int k=0; k<numFacBindingStates; k++){
							MUs[j][m] = k==0 && i==0? Qijk[i][j][k]*Xs[i][m]: MUs[j][m]+Qijk[i][j][k]*Xs[i][m];
							denMUs[j][m] = denMUs[j][m]+Qijk[i][j][k];
						}
					}
				}
			}
			
			//Normalize
			for(int j=0; j<numChromStates; j++){
				for(int m=0; m<M; m++){
					MUs[j][m] = MUs[j][m]/denMUs[j][m];
				}
			}
		}
		

		//--------------------SIGMAc update --------------------------------------
		
		//Compute
		double[][] denSIGMAc = new double[numChromStates][C];
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				for(int i=0; i<N; i++){
					for(int k=0; k<numFacBindingStates; k++){
						SIGMAc[j][c] = k==0 && i==0? Qijk[i][j][k]*Math.pow(((double)Xc[i][c] - MUc[j][c]) , 2.0): SIGMAc[j][c]+Qijk[i][j][k]*Math.pow(((double)Xc[i][c] - MUc[j][c]) , 2.0);
						denSIGMAc[j][c] = denSIGMAc[j][c]+Qijk[i][j][k];
					}
				}
			}
		}
		
		//Normalize and taking the square root
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				SIGMAc[j][c] = Math.sqrt(SIGMAc[j][c]/denSIGMAc[j][c]);
			}
		}
		
		//------------------SIGMAf update------------------------------------------
		
		//Compute
		double[][] denSIGMAf = new double[numFacBindingStates][F];
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				for(int i=0; i<N; i++){
					for(int j=0; j< numChromStates; j++){
						SIGMAf[k][f] = j==0 && i==0? Qijk[i][j][k]*Math.pow(((double)Xf[i][f] - MUf[k][f]) , 2.0) : SIGMAf[k][f]+Qijk[i][j][k]*Math.pow(((double)Xf[i][f] - MUf[k][f]) , 2.0);
						denSIGMAf[k][f] = denSIGMAf[k][f]+Qijk[i][j][k];
					}
				}
			}
		}
		
		//Normalize and taking the square root
		for(int k=0; k<numFacBindingStates; k++){
			for(int f=0; f<F; f++){
				SIGMAf[k][f] = Math.sqrt(SIGMAf[k][f]/denSIGMAf[k][f]);
			}
		}
		
		
		//--------------------SIGMAs update, if in seqState ------------------------
		
		//Compute
		if(this.seqState){
			double[][] denSIGMAs = new double[numChromStates][M];
			for(int j=0; j<numChromStates; j++){
				for(int m=0; m<M; m++){
					for(int i=0; i<N; i++){
						for(int k=0; k<numFacBindingStates; k++){
							SIGMAs[j][m] = k==0 && i==0? Qijk[i][j][k]*Math.pow(((double)Xs[i][m] - MUs[j][m]) , 2.0): SIGMAs[j][m]+Qijk[i][j][k]*Math.pow(((double)Xs[i][m] - MUs[j][m]) , 2.0);
							denSIGMAs[j][m] = denSIGMAs[j][m]+Qijk[i][j][k];
						}
					}
				}
			}
			
			//Normalize and taking the square root
			for(int j=0; j<numChromStates; j++){
				for(int m=0; m<M; m++){
					SIGMAs[j][m] = Math.sqrt(SIGMAs[j][m]/denSIGMAs[j][m]);
				}
			}
		}
		
		//---------------------Bjk update -------------------------------------------
		
		//Compute
		double[] denBjk = new double[numFacBindingStates];
		for(int k=0; k<numFacBindingStates; k++){
			for(int j=0; j< numChromStates; j++){
				for(int i=0; i<N; i++){
					Bjk[j][k] = i==0 ? Qijk[i][j][k] : Bjk[j][k]+Qijk[i][j][k];
					denBjk[k] = denBjk[k]+Qijk[i][j][k];
				}
			}
		}
		
		//Normalize
		for(int k=0; k<numFacBindingStates; k++){
			for(int j=0; j< numChromStates; j++){
				Bjk[j][k] = Bjk[j][k]/denBjk[k];
			}
		}
	}
	
	
	//setters
	public void setSeqMode(Sequences seqs, double[][] Xs){
		this.seqState = true;
		this.setXs(Xs);
		this.M = Xs[0].length;
		this.setSequences(seqs);
		this.setInitialSeqParams();
		this.trainMUs = new double[this.total_itrs+1][this.numChromStates][this.M];
		this.trainSIGMAs = new double[this.total_itrs+1][this.numChromStates][this.M];
		
	}
	private void setSequences(Sequences seqs){this.seqs = seqs;}
	private void setXs(double[][] Xs){this.Xs = Xs;}
	private void setInitialSeqParams(){this.initializeSeqParams();}
	public void serTotalNoItrs(int n){
		if(this.total_itrs == 0){
			this.total_itrs = n;
		}
	}
	
	//Accessors
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
