package org.seqcode.projects.shaun;

import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.probability.NormalDistribution;

/**
 * Make a Gaussian distribution for GPS
 * @author mahony
 *
 */
public class MakeGaussianDistrib {

	public static void main(String[] args){
		int start=-100, end=100;
		Double mean = 0.0;
		Double sigma = 5.0;
		if(Args.parseArgs(args).contains("mean")){
			mean = Args.parseDouble(args, "mean", 0);
		}
		if(Args.parseArgs(args).contains("sigma")){
			sigma = Args.parseDouble(args, "sigma", 5);
		}
		
		NormalDistribution gaussian = new NormalDistribution(mean, sigma*sigma);
		for(int i=start; i<=end; i++){
			System.out.println(i+"\t"+gaussian.calcProbability((double)i));
		}
	}
}
