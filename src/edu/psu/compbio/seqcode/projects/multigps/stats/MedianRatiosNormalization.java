package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.util.ArrayList;
import java.util.Collections;

/**
 * MedianRatiosNormalization: a Normalization class that implements the median of ratios normalization
 * as proposed by Anders & Huber, Genome Biology, 2010
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class MedianRatiosNormalization extends Normalization{

	int ref=0;
	
	/**
	 * Constructor: provide the number of samples, the fraction of extreme M values to trim, and the fraction of extreme A values to trim
	 * @param numSamples
	 * @param Mtrim
	 * @param Atrim
	 */
	public MedianRatiosNormalization(int numSamples) {
		super(numSamples);
	}

	@Override
	public double[] normalize(CountsDataset data) {
		
		//Set one sample as the reference (the deepest sequenced sample in the focal condition)
		double maxTotal=0;
		for(int s=0; s<data.numSamples; s++)
			if(data.design[s] == data.focalCondition)
				if(data.totals[s]>maxTotal){
					ref=s; maxTotal=data.totals[s];
				}
		
		//Make a pseudo-reference
		double [] pseudo = new double[data.numUnits];
		double frac = 1/(double)data.numSamples;
		for(int d=0; d<data.numUnits; d++){
			double prod = 1;
			for(int s=0; s<data.numSamples; s++)
				prod*=data.getCount(d,s);
			if(prod==0)
				pseudo[d]=-1;
			else
				pseudo[d] = Math.pow(prod, frac);
		}
		//Arrays of ratios against the pseudo reference
		for(int s=0; s<data.numSamples; s++){
			ArrayList<Double> ratios = new ArrayList<Double>(); 
			for(int d=0; d<data.numUnits; d++)
				if(pseudo[d]>0)
					ratios.add(data.getCount(d,s)/pseudo[d]);
			//Scaling factor is the median ratio
			Collections.sort(ratios);
			
			depthScaling[s]= (data.totals[s] / data.totals[ref]);
			scalingFactors[s] = ratios.get(ratios.size()/2);
			propScaling[s] = scalingFactors[s]/depthScaling[s];
		}
		data.setScalingFactors(scalingFactors);
		return scalingFactors;
	}

	@Override
	public String selfDescriptor() {
		return "Median ratios normalization as proposed by Anders & Huber, Genome Biology, 2010 (DESeq).";
	}

}
