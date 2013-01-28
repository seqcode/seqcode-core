package edu.psu.compbio.seqcode.projects.multigps.framework;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * PoissonBackgroundModel: A BackgroundModel using the Poisson model
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class PoissonBackgroundModel extends BackgroundModel{

	double lambda; 
	
	public PoissonBackgroundModel(int modelType, double logconfidence, double totalReads, double genomeLength, double mappableGenome, double binWidth, char strand, double scaling, boolean useThisExpt) {
		super(modelType, logconfidence, totalReads, genomeLength, mappableGenome, binWidth, strand, scaling, useThisExpt);
		countThreshold = calcCountThreshold();
	}

	//Does the hit count in a region pass the threshold? 
	public boolean passesThreshold(int count) {
		if(count>=countThreshold)
			return true;
		else
			return false;
	}
	//Does the hit count in a region stay under the threshold? 
	public boolean underThreshold(int count) {
		if(count<countThreshold)
			return true;
		else
			return false;
	}

	//Set Poisson thresholds
	protected int calcCountThreshold(){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		lambda = (totalReads*binWidth)/(regionLength*mappableRegion); 
		P.setMean(lambda);
		double l=1;
		for(int b=1; l>confThreshold; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		return(Math.max(1,countThres));
	}
}

