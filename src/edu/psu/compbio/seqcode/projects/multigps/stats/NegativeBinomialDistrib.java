package edu.psu.compbio.seqcode.projects.multigps.stats;

import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;
import cern.jet.random.Gamma;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * NegativeBinomial: distribution of the number of failures (X) before the rth success in independent trials, with success 
 * probability p in each trial. 
 * 
 * This class implements a random number generator that samples from the NB. PDF & CDF will be added/ 
 * 
 * This class is necessary because of a bug in the Colt implementation of NegativeBinomial (in the nextInt method), and also because
 * the Colt implementation does not allow a real-valued r. 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class NegativeBinomialDistrib {
	protected double r;	//Number of successes 
	protected double p;	//Probability of success
	protected Gamma gamma; //Gamma distribution
	protected Poisson poisson; //Poisson
	protected NegativeBinomialDist nb; //SSJ NegativeBinomial
	
	/**
	 * NegativeBinomial (r p parameterization)
	 * @param r : Number of successes 
	 * @param p : Probability of success in each trial
	 */
	public NegativeBinomialDistrib(double r, double p){
		this.r = r;
		this.p = p;
		gamma = new Gamma(1, 1.0, new DRand());
		poisson = new Poisson(10, new DRand());
		nb = new NegativeBinomialDist(r, p);
	}
	
	/**
	 * Update the r & p parameters direcctly
	 * @param r
	 * @param p
	 */
	public void setRandP(double r, double p){
		this.r = r;
		this.p = p;
	}
	
	/**
	 * Update the r & p parameters via the mean and variance 
	 * @param mean
	 * @param var
	 */
	public void setMeanVar(double mean, double var){
		r = (mean*mean)/(var-mean);
		p = mean/var;
	}
	
	/**
	 * Returns a random number from the distribution.
	 * @return integer
	 */
	public int nextInt(){
		return nextInt(r, p);
	}

	/**
	 * Returns a random number from the distribution, bypassing internal state
	 * This method samples a random number from the NegativeBinomial
	 * distribution with parameters r (no. of successes in n independent trials) 
	 * and p (probability of success)
	 * Valid for r > 0, 0 < p < 1.
	 * If G from Gamma(r) then K from Poiss((1-p)G/p) is NB(r,p)--distributed. 
	 * REFERENCE: - J.H. Ahrens, U. Dieter (1974): Computer methods for
	 * sampling from gamma, beta, Poisson and * binomial distributions,
	 * Computing 12, 223--246.
	 * 
	 * @return integer
	 */
	public int nextInt(double nb_r, double nb_p){
		double x = (1.0 - nb_p)/nb_p;
		double y = x * gamma.nextDouble(nb_r, 1.0);
		return poisson.nextInt(y);
	}
	
	/**
	 * Call SSJ NegativeBinomial pdf
	 * @param x Number of failures before the r-th success
	 * @return
	 */
	public double pdf(int x){
		return nb.prob(x);
	}
	
	/**
	 * Call SSJ NegativeBinomial static pdf
	 * @param x Number of failures before the r-th success
	 * @param nb_r Number of successes
	 * @param nb_p Probability of success
	 * @return
	 */
	public static double pdf(int x, double nb_r, double nb_p){
		return NegativeBinomialDist.prob(nb_r, nb_p, x);
	}
	
	/**
	 * Call SSJ NegativeBinomial cdf
	 * @param x Number of failures before the r-th success
	 * @return
	 */
	public double cdf(int x){
		return nb.cdf(x);
	}
	
	/**
	 * Call SSJ NegativeBinomial static cdf
	 * @param x Number of failures before the r-th success
	 * @param nb_r Number of successes
	 * @param nb_p Probability of success
	 * @return
	 */
	public static double cdf(int x, double nb_r, double nb_p){
		return NegativeBinomialDist.cdf(nb_r, nb_p, x);
	}
}
