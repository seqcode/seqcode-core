package org.seqcode.ml.classification;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.trees.*;

public class WekaRandomForest {
	
	public static void main(String[] args){
		RandomForest rf = new RandomForest();
		AbstractClassifier.runClassifier(rf, args);
	}
}
