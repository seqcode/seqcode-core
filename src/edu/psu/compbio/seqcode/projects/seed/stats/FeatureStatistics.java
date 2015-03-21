package edu.psu.compbio.seqcode.projects.seed.stats;

import java.util.Collections;
import java.util.List;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.projects.seed.features.Feature;

public class FeatureStatistics {

	protected Binomial binomial;
	
	public FeatureStatistics(){
		binomial = new Binomial(100,.5, new DRand());
	}
	
	
	/**
	 * Binomial CDF assuming scaled control. Uses COLT binomial test.
	 * Tests equality of signal & scaled control counts. 
	 * @param k = scaled control
	 * @param n = scaled control+signal
	 * @param minFoldChange = minimun fold difference between signal & control
	 * @return
	 */
	public double binomialPValue(double k, double n, double minFoldChange){
		double pval=1;
		synchronized(binomial){
			binomial.setNandP((int)Math.ceil(n), 1.0 / (minFoldChange + 1.0));
			pval = binomial.cdf((int) Math.ceil(k));
		}
        return(pval);		
	}
	
	/**
	 * In-place FDR-based multiple testing correction by computing q-values. 
	 *  Benjamini-Hochberg procedure, enforcing monotonicity of p-value increases as defined in: 
	 *  Storey JD. A direct approach to false discovery rates. Journal of the Royal Statistical Society. 2002;64:479-498
	 * 
	 */
	public void benjaminiHochbergCorrection(List<Feature> features){
		if(features.size()>0 && features.get(0).scoreIsAPValue){
			double total = (double)features.size();
			double rank = (double)features.size();
			
			//Correct the lowest ranked item first
			double q = Math.min(features.get(features.size()-1).getScore()*(total/rank), 1);
			features.get(features.size()-1).setScore(q);
			//Iterate backwards
			for(int i=features.size()-2; i>=0; i--){
				Feature fi = features.get(i);
				rank--;
				//Set q-value based 
				q = Math.min(fi.getScore()*(total/rank), features.get(i+1).getScore());
				fi.setScore(q);
			}
		}
	}
	
	
	//TODO: IDR method
	
}
