package edu.psu.compbio.seqcode.projects.naomi.multiscalesignal;

import java.util.Arrays;

/**
 * Hierarchical Hidden Markov Model (HHMM)
 *
 * Kevin Murphy (2001) Hierarchical HMMs
 * Shai Fine et al (1998) The Hierarchical Hidden Markov Model: Analysis and Applications
 * 
 * @author naomi yamada
 *
 **/

public class HierarchicalHMM {

	protected int numStates; 
	protected int numObs;
	protected int depth; //hierarchy depth
	protected double o [][]; //observations	o[numObs][depth]
	protected int q [][]; //ith state at level d ; temporary set as int to do the initial development
	protected double pi [][]; // initial state distributions
	protected double a [][][]; // state transition probabilities
	protected double b [][][]; //emission probabilities	
	
	/**
	 * Class constructor. Initializes a HHMM given a number of states and depth.
	 * 
	 * @param numStates		number of states
     * @param depth			depth
	 **/	
	public HierarchicalHMM (int numStates, int depth , int numObs){
		
		if (numStates <1) numStates = 1;
		if (depth <1) depth = 1;
		this.numStates = numStates;
		this.numObs = numObs;
		this.depth = depth;
		
		pi = new double [numStates][depth];
		a = new double [numStates][numStates][depth];
		b = new double [numStates][numObs][depth];
		
	}

	/**
	 * Baum-Welch algorithm for hidden Markov models. Given an observation of mean signal,
	 * it will train HHMM using output to increase probability of this HHMM generating output.
	 * 
	 * @param o		signal mean
	 * @param itr	the number of iterations performed
	 **/	
	public void baumWelch (int[][] o, int itr){
		
		int T = o.length;
		double [][][] forward;
		double [][][] backward;
		double [][][] downward;
		double [][][] upward;
		
		double pi1[] = new double [numStates];
		double a1[][][] = new double [numStates][numStates][depth];
		double b1[][][] = new double [numStates][numObs][depth];
		
		for (int s = 0; s <itr; s++){
			//calculate forward, backward, downward, upward variables from the current model
//			forward = forwardProc(o);
//			backward = backwardProc(o);
//			upward = upwardProc(o);
//			downward = downwardProc(o);
			
			//re-estimation of initial state probabilities
			for (int i = 0; i <numStates; i++){
	//			pi1[i] = gamma(i,0,o,forward,backward);
			}
			
			//re-estimation of transition probabilities
			
		      // re-estimation of transition probabilities
		      for (int i = 0; i < numStates; i++) {
				for (int j = 0; j < numStates; j++) {
				  double num = 0;
				  double denom = 0;
				  for (int t = 0; t <= T - 1; t++) {
//				    num += xi(t, i, j, o, forward, backward);
	//			    denom += gamma(i, t, o, forward, backward);
				  }
	//			  a1[i][j] = divide(num, denom);
				}
		      }
			
			//re-estimation of emission probabilities
			
			
		}
		
	}

	/**
	 * Calculation of forward variables f(i,t,d) for state i at time t at depth d
	 * for mean o with the current HHMM parameters.
	 * 
	 * @param o		signal mean
	 * @return 		a 3d array containing forward variables f(i,t,d) over states, time, and depth.
	 **/			
	public double [][][] forwardProc(double[][] o){
		
		int T = o[0].length;
		double [][][] forward = new double [numStates][T][depth];
		
		//Bottom level: d=D
		//for each state, calculate probability based on state initial distributions of d-1
		
		for (int i = 0; i <numStates; i++)
			forward[i][0][depth] = pi[q[i][depth-1]][depth-1]; //multiplied by observation
		
		
		
		// initialization (time 0)
		for (int i = 0; i <numStates; i++)
			forward[i][0][0] = pi[i][0] *b[i][0][0];
		
		// induction when depth is zero
		for (int t = 0 ; t <= T-2 ; t++){
			for (int j = 0 ; j < numStates ; j++){
				forward[j][t+1][0] = 0;
				for (int i = 0; i< numStates; i++){
					forward[j][t+1][0] += (forward[i][t][0]*a[i][j][0]);
				forward[j][t+1][0] *= b[j][t+1][0];
				}
			}
		}

		
		return forward;
		
	}
	
	public double [][][] backwardProc(double[][] o){
		
		int T = o[0].length;
		double [][][] backward = new double [numStates][T][depth];
		
		// initialization (time 0)
		for (int i = 0; i <numStates; i++)
			backward[i][T-1][0] = 1;
		
		// induction at depth zero
		for (int t = T -2 ; t >= 0; t--){
			for (int j = 0; j < numStates ; j ++){
				backward[j][t][0] = 0;
				for (int i = 0; i < numStates; i++){
					backward[j][t][0] += (backward[j][t+1][0] * a[i][j][0] * b[j][t+1][0]);
				}
			}
		}
		
		/**
		//induction
		for (int d = depth ; d >= 0; d--){		
			for (int t = T -2 ; t >= 0; t--){
				for (int j = 0; j < numStates ; j ++){
					backward[j][t][d] = 0;
					for (int i = 0; i < numStates; i++){
						backward[j][t][d] += (backward[j][t+1][d] * a[i][j][d] * b[j][t+1][d]);
					}
				}
			}
		}
		**/
		
		return backward;
	}
	
	public double [][][] upwardProc(double[][] o){
		
		int T = o[0].length;
		double [][][] upward = new double [numStates][T][depth];
		
		return upward;
	}
	
	public double [][][] downwardProc(double[][] o){
		
		int T = o[0].length;
		double [][][] downward = new double [numStates][T][depth];
		
		return downward;
	}

	/**
	   * Calculation of xi_t(i, j), which is the probability
	   * P(i_t = s_i, i_t+1 = s_j | o, hmm), that is, 
	   * the probability of being in state i at time t and state j 
	   * at time t+1 given observation sequence o and this HMM.
	   *
	   * @param t		the time
	   * @param i		the number of state s_i
	   * @param j		the number of state s_j
	   * @param o		the observation sequence
	   * @param fwd		the Forward variables for o
	   * @param bwd		the Backward variables for o
	   * @return		P(i_t = s_i, i_t+1 = s_j | o, hmm)
	   */
	public double xi(int t, int i, int j, double[][] o, double[][][] forward, double[][][] backward){
	    
		 double num, denom = 0.0;
		 
		 // numerator
		 if (t == o.length - 1)
			 num = forward[i][t][0] * a[i][j][0];
		 else
			 num = forward[i][t][0] * a[i][j][0] * b[j][t+1][0] * backward[j][t+1][0];

		 // denominator
		 for (int k = 0; k < numStates; k++)
			 denom += (forward[k][t][0] * backward[k][t][0]);

		 return divide(num, denom);
		 }	
	


	  /**
	   * Calculation of gamma_t(i), which is the probability
	   * P(i_t = s_i | o, hmm), that is, the probability
	   * of being in state i at time given observation sequence 
	   * o and this HMM.
	   *
	   * @param i		the number of state s_i
	   * @param t		the time
	   * @param o		the observation sequence
	   * @param fwd		the Forward variables for o
	   * @param bwd		the Backward variables for o
	   * @return		P(i_t = s_i | o, hmm)
	   */
	public double gamma_in(int i, int t, double[][] o, double[][][] forward, double[][][] backward){
		  
	    double num, denom = 0.0;
	    
	    // numerator
	    num = forward[i][t][0] * backward[i][t][0];

		// denominator
	    for (int j = 0; j < numStates; j++)
	      denom += forward[j][t][0] * backward[j][t][0];

	    return divide(num, denom);
	  }		
	
	public double gamma_out(int i, int t, double[][] o, double[][][] forward, double[][][] backward){
		  
	    double num, denom = 0.0;
	    
	    // numerator
	    num = forward[i][t][0] * backward[i][t][0];

		// denominator
	    for (int j = 0; j < numStates; j++)
	      denom += forward[j][t][0] * backward[j][t][0];

	    return divide(num, denom);
	  }		
	
	public double eta_in(int i, int t, double[][] o, double[][][] forward, double[][][] backward){
		  
	    double num, denom = 0.0;
	    
	    // numerator
	    num = forward[i][t][0] * backward[i][t][0];

		// denominator
	    for (int j = 0; j < numStates; j++)
	      denom += forward[j][t][0] * backward[j][t][0];

	    return divide(num, denom);
	  }	
	
	public double eta_out(int i, int t, double[][] o, double[][][] forward, double[][][] backward){
		  
	    double num, denom = 0.0;
	    
	    // numerator
	    num = forward[i][t][0] * backward[i][t][0];

		// denominator
	    for (int j = 0; j < numStates; j++)
	      denom += forward[j][t][0] * backward[j][t][0];

	    return divide(num, denom);
	  }	

	
	public double chi(int t, int i, double[][] pi, double[][][] forward, double[][][] backward){
		
		double num, denom = 0.0;
		
		// numerator
		if (t == o.length - 1)
			num = forward[i][t][0] * a[i][t][0];
		else
			num = forward[i][t][0] * a[i][t][0] * b[i][t+1][0] * backward[i][t+1][0];

		// denominator
		for (int k = 0; k < numStates; k++)
			denom += (forward[k][t][0] * backward[k][t][0]);
		
		return divide(num, denom);
	}
		
	
	public int findIndex(double [][] array, int d, double value){
		for (int i = 0; i < array[0].length; i++){
			if (array[i][0] == value)
				return i;
		}
		return -1;
	}
	
	/* Divides two doubles.
	 * 0 / 0 = 0!
	 */
	private static double divide(double n, double d){
		if (n == 0)
			return 0;
		else
			return n / d;
	  }
	
	
}
