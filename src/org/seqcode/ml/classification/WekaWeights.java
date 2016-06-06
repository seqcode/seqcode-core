package org.seqcode.ml.classification;

import weka.classifiers.*;

public class WekaWeights {
	
	public static void main(String[] args) throws Exception{
		if(args.length != 1){
			System.err.println("give model file as param");
			System.exit(1);
		}
		
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(args[0]);
		System.err.println("Reading "+args[0]+":"+cls.getClass().getName());
		String[] attributes = null;
		double[] weights = null;
		if(cls instanceof BaggedRandomForest){
			weights = ((BaggedRandomForest) cls).getAttributeWeights();
			attributes = ((BaggedRandomForest) cls).getAttributes();
		}else
			throw new Error("unknown classifier");
		
		for(int i=0; i<weights.length; i++){
			System.out.println(attributes[i]+" "+weights[i]);
		}
		
	}

}
