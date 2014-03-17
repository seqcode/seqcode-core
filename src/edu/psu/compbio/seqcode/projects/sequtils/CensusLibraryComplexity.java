package edu.psu.compbio.seqcode.projects.sequtils;

import jsc.distributions.LogarithmicSeries;
import jsc.distributions.Lognormal;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariateOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.multigps.stats.NegativeBinomialDistrib;

/**
 * Census method for estimating library complexity.
 * 
 * Adapted from Matt Edwards (https://github.com/matted/census).
 * 
 * Fits a truncated Poisson and Negative Binomial to a histogram of observed read counts per base.
 * 
 * @author mahony
 *
 */
public class CensusLibraryComplexity {

	RealValuedHistogram histo; //Histogram of per-base read counts
	double xsum=0, xcount=0;
	int truncLowerBound = 1;
	int truncUpperBound = 100; 
	double relativeAccuracy = 1.0e-6;
	double absoluteAccuracy = 1.0e-4;
	boolean verbose=false;
	double estimatedPoissonLambda=0;
	double [] estimatedNegBinomialRP={0,0};
	double [] estimatedLogNormMS={0,0};
	double estimatedPoissonLibSize=0;
	double estimatedNegBinomialLibSize=0;
	double estimatedLogNormLibSize=0;
	
	public CensusLibraryComplexity(RealValuedHistogram h, int truncLower, int truncUpper){
		histo = h;
		truncLowerBound = truncLower;
		truncUpperBound = truncUpper;
	}
	
	//Accessors
	public void setVerbose(boolean v){verbose=v;}
	public double getEstimatedPoissonLambda(){return estimatedPoissonLambda;}
	public double[] getEstimatedNegBinomialRP(){return estimatedNegBinomialRP;}
	public double[] getEstimatedLogNormMS(){return estimatedLogNormMS;}
	public double getEstimatedPoissonLibSize(){return estimatedPoissonLibSize;}
	public double getEstimatedNegBinomialLibSize(){return estimatedNegBinomialLibSize;}
	public double getEstimatedLogNormLibSize(){return estimatedLogNormLibSize;}
	public double getEstimatedNegBinomialMean(){ //r*(1-p)/p;
		return estimatedNegBinomialRP[0] * ((1.0-estimatedNegBinomialRP[1])/ estimatedNegBinomialRP[1]);
	}
	public double getEstimatedNegBinomialVariance(){ //mean + k*mean^2
		double mean = getEstimatedNegBinomialMean();
		double k = getEstimatedNegBinomialGammaK();
		return (mean + k*mean*mean);
	}
	//alpha parameter from gamma distribution formulation
	public double getEstimatedNegBinomialAlpha(){return estimatedNegBinomialRP[0];}
	//beta parameter from gamma distribution formulation
	public double getEstimatedNegBinomialBeta(){return (1.0-estimatedNegBinomialRP[1])/ estimatedNegBinomialRP[1];}
	//squared coefficient of variation of the gamma distribution that describes the rates ("sequencing probabilities") of the molecules in the library
	public double getEstimatedNegBinomialGammaK(){return 1/estimatedNegBinomialRP[0];} 
	//coefficient of variation of the gamma distribution that describes the rates ("sequencing probabilities") of the molecules in the library
	public double getEstimatedNegBinomialGammaCV(){return 1/(Math.sqrt(estimatedNegBinomialRP[0]));}
	//Mean of log normal
	public double getEstimatedLogNormMean(){ return new Lognormal(estimatedLogNormMS[0], estimatedLogNormMS[1]).mean();}
	//Coverage of distinct observations
	public double getEstimatedPoissonObservedCoverage(){
		DRand re = new DRand();
		Poisson poiss = new Poisson(estimatedPoissonLambda, re);
		return (1.0 - poiss.pdf(0));
	}
	//Coverage of distinct observations
	public double getEstimatedNegBinomialObservedCoverage(){
		NegativeBinomialDistrib nb = new NegativeBinomialDistrib(estimatedNegBinomialRP[0], estimatedNegBinomialRP[1]);
		return (1.0 - nb.pdf(0));
	}
	//Coverage of distinct observations
	public double getEstimatedLogNormObservedCoverage(){
		Lognormal ln = new Lognormal(estimatedLogNormMS[0], estimatedLogNormMS[1]);
		return (1.0 - ln.cdf(1));
	}
	//number of distinct observations given a coverage count
	public double getEstimatedPoissonObservedGivenCount(double readCount){
		DRand re = new DRand();
		Poisson poiss = new Poisson(readCount/estimatedPoissonLibSize, re);
		return (1.0 - poiss.pdf(0));
	}
	//number of distinct observations given a coverage count
	public double getEstimatedNegBinomialObservedGivenCount(double readCount){
		//the model for future prediction is that alpha (i.e. r) is constant since it's a coverage-independent property of the library, and you increase depth by cranking up the Poisson mean.
		double new_mean = readCount/estimatedNegBinomialLibSize;
		double old_beta = (1.0-estimatedNegBinomialRP[1])/ estimatedNegBinomialRP[1];
		double new_beta = old_beta * new_mean / (xsum / estimatedNegBinomialLibSize); //scale beta to stay proportional to new mean
		NegativeBinomialDistrib nb = new NegativeBinomialDistrib(estimatedNegBinomialRP[0], 1.0 - new_beta/(1.0+new_beta));
		return (1.0 - nb.pdf(0));
	}
	//number of distinct observations given a coverage count
	public double getEstimatedLogNormObservedGivenCount(double readCount){
		
		return(0);
	}
	
	//Execute the estimation methods
	public void execute(){
		DRand re = new DRand();
		int left=truncLowerBound, right=truncUpperBound;
		xsum=0; xcount=0;
		for(double i=left; i<=right; i++){
			xsum += i*histo.getBin( histo.getBinContainingVal(i));
			xcount += histo.getBin( histo.getBinContainingVal(i));
		}
		double xavg = xsum/xcount;
		
		//Fit a Poisson
		UnivariateFunction func = new truncPoisson(xavg, left, right);
		UnivariateOptimizer solver = new BrentOptimizer(relativeAccuracy, absoluteAccuracy);
		UnivariatePointValuePair pvp = solver.optimize(new MaxEval(100), new UnivariateObjectiveFunction(func), GoalType.MINIMIZE, new SearchInterval(0.001, 200.0));
		double lambda = pvp.getPoint();
		estimatedPoissonLambda=lambda;
		//if(verbose)
		//	System.out.println(String.format("Poisson: xavg= %.5f\tlambda= %.5f", xavg, lambda));
		
		//Fit a negative binomial
		//Iterate over several initialization choices to find the best one (or at least a better one):
		double bestNegLL = Double.POSITIVE_INFINITY;
		double [] bestRP={0,0};
		double [] ks = {0.01, 0.1, 0.5, 0.9};
		truncNegativeBinomial truncNB = new truncNegativeBinomial(histo, left, right);
		for(int i=0; i<ks.length; i++){
			double mean = lambda; 
			double k = ks[i];
			double var = mean + k*mean*mean;
			double r = (mean*mean)/(var-mean); 
			double p = mean/var;
			double [] initRP = {r, p};
			double[] lb = {r/10, p/10}; double[] up = {r*10, 1.0-1e-9};
			
			//if(verbose)
			//	System.out.println(String.format("initRP: %.5f %.5f",  initRP[0], initRP[1]));
			MultivariateFunction mfunc = truncNB;
			MultivariateOptimizer msolver = new BOBYQAOptimizer(5);
			try{
				PointValuePair mpvp = msolver.optimize(new MaxEval(10000), new ObjectiveFunction(mfunc), GoalType.MINIMIZE, new InitialGuess(initRP), new SimpleBounds(lb, up) );
				double[] currRP = mpvp.getPoint();
				double currNegLL = truncNB.value(currRP);
				if(currNegLL < bestNegLL){
					bestNegLL= currNegLL;
					bestRP = currRP;
				}
			}catch(TooManyEvaluationsException e){
				
			}catch(MathIllegalStateException e){
				
			}
		}
		estimatedNegBinomialRP=bestRP;
		//if(verbose)
		//	System.out.println(String.format("NegativeBinomial: xavg= %.5f\tR= %.5f\tP= %.5f\tmean=%.5f", xavg, bestRP[0], bestRP[1], bestRP[0]*(1-bestRP[1])/bestRP[1]));
		
		//Fit a log normal
		bestNegLL = Double.POSITIVE_INFINITY;
		double [] bestMS={0,0};
		truncLogNormal truncLN = new truncLogNormal(histo, left, right);
		MultivariateFunction mfunc = truncLN;
		MultivariateOptimizer msolver = new BOBYQAOptimizer(5);
		double[] lb = {0, 1e-6}; double[] up = {10, 10};
		try{
			PointValuePair mpvp = msolver.optimize(new MaxEval(10000), new ObjectiveFunction(mfunc), GoalType.MINIMIZE, new InitialGuess(new double[]{0.0, 0.5}), new SimpleBounds(lb, up) );
			double[] currMS = mpvp.getPoint();
			double currNegLL = truncLN.value(currMS);
			if(currNegLL < bestNegLL){
				bestNegLL= currNegLL;
				bestMS = currMS;
			}
		}catch(TooManyEvaluationsException e){
			
		}
		estimatedLogNormMS = bestMS;
		
		//Estimate library sizes
		Poisson poiss = new Poisson(lambda, re);
		estimatedPoissonLibSize = (xcount / (poiss.cdf(right) - poiss.cdf(left - 1)));
		NegativeBinomialDistrib negBinom = new NegativeBinomialDistrib(bestRP[0], bestRP[1]);
		estimatedNegBinomialLibSize = (xcount / Math.max(1e-12, negBinom.cdf(right) - negBinom.cdf(left - 1)));
		Lognormal lognor = new Lognormal(estimatedLogNormMS[0], estimatedLogNormMS[1]);
		estimatedLogNormLibSize = (xcount / lognor.cdf(right) - lognor.cdf(left - 1));
		
		//Print observed and estimated histos
		if(verbose){
			System.out.println("\nObserved/Expected per-base counts:\nCount\tObs\tPoissonExp\tNegBinomialExp\tLogNormal");
			for(int i=0; i<histo.getHistoStop(); i++){
				double obs = histo.getBin(histo.getBinContainingVal(i));
				double poissExp = estimatedPoissonLibSize * poiss.pdf(i);
				double nbExp = estimatedNegBinomialLibSize * negBinom.pdf(i);
				double lnExp = estimatedLogNormLibSize * (lognor.cdf(i+1)-lognor.cdf(i));
				System.out.println(String.format("%d\t%.0f\t%.0f\t%.0f\t%.0f",i, obs, poissExp, nbExp,lnExp));
				//System.out.println(String.format("%d\t%.0f\t%.0f\t%.0f",i, obs, poissExp, nbExp));
			}
		}
	}
	protected double logTrunc(double x){return Math.log(Math.max(x, 1e-300));}
	
	
	//A truncated Poisson function
	private class truncPoisson implements UnivariateFunction {
		protected DRand re = new DRand();
		protected double xavg;
		protected int left, right;
		
		public truncPoisson(double xavg, int left, int right){
			this.xavg = xavg;
			this.left = left;
			this.right = right;
		}
		public double value(double L){
			//if(verbose)
		    //    System.out.println(String.format("val: %.3f %.5f %.5f",L, -(-logTrunc(K(L, left, right)) - L + xavg * logTrunc(L)), -logTrunc(K(L, left, right))));
		    return -(-logTrunc(K(L, left, right)) - L + xavg * logTrunc(L));		
		}
		public double K(double L, int left, int right){
			Poisson poiss = new Poisson(L, re); 
		    return poiss.cdf(right) - poiss.cdf(left - 1);
		}
	}
	
	//A truncated Negative Binomial function
	private class truncNegativeBinomial implements MultivariateFunction {
		protected RealValuedHistogram hist;
		protected int left, right;
		
		public truncNegativeBinomial(RealValuedHistogram h, int left, int right){
			this.hist = h;
			this.left = left;
			this.right = right;
		}
		public double value(double[] rp){
			double LL = 0;
			
			double logK = logTrunc(NegativeBinomialDistrib.cdf(right, rp[0], rp[1]) - NegativeBinomialDistrib.cdf(left-1, rp[0], rp[1]));
    	    for(int index=left; index<=right; index++){
    	        if(index >=hist.getHistoStart() && index <=hist.getHistoStop())
    	        	LL += (logTrunc(NegativeBinomialDistrib.pdf(index,rp[0], rp[1])) - logK) * histo.getBin(histo.getBinContainingVal(index));
    	    }
    	    //if(verbose)
			//	System.out.println(String.format("val: %.5f %.5f %.5f %.5f",  rp[0], rp[1], -LL, -logK));
    	    
    	    return LL==0.0 ? Double.MAX_VALUE : -LL;
		}
		public double K(double r, double p, int left, int right){
			NegativeBinomialDistrib negbinom = new NegativeBinomialDistrib(r, p); 
		    //if(verbose)
		    //    System.out.println(String.format("K: %.5f %.5f %.5f %.5f %d %d",negbinom.cdf(right),negbinom.cdf(left - 1), n, p, left, right));
		    return negbinom.cdf(right) - negbinom.cdf(left - 1);
		}
	}

	//A truncated Log-series function
	private class truncLogNormal implements MultivariateFunction {
		protected RealValuedHistogram hist;
		protected int left, right;
		
		public truncLogNormal(RealValuedHistogram h, int left, int right){
			this.hist = h;
			this.left = left;
			this.right = right;
		}
		public double value(double[] ms){
			double LL = 0;
			Lognormal lognor = new Lognormal(ms[0], ms[1]);
			
			double logK = logTrunc(lognor.cdf(right) - lognor.cdf(left-1));
    	    for(int index=left; index<=right; index++){
    	        if(index >=hist.getHistoStart() && index <=hist.getHistoStop())
    	        	LL += (logTrunc(lognor.pdf(index)) - logK) * histo.getBin(histo.getBinContainingVal(index));
    	    }
    	    //if(verbose)
			//	System.out.println(String.format("val: %.5f %.5f %.5f %.5f",  ms[0], ms[1], -LL, -logK));
    	    
    	    return LL==0.0 ? Double.MAX_VALUE : -LL;
		}
		public double K(double m, double s, int left, int right){
			Lognormal lognor = new Lognormal(m, s); 
		    //if(verbose)
		    //    System.out.println(String.format("K: %.5f %.5f %.5f %.5f %d %d",negbinom.cdf(right),negbinom.cdf(left - 1), n, p, left, right));
		    return lognor.cdf(right) - lognor.cdf(left - 1);
		}
	}
}
