package org.seqcode.ml.classification;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.bayes.*;

public class WekaNaiveBayes {
	
	public static void main(String[] args){
		NaiveBayes nb = new NaiveBayes();
		AbstractClassifier.runClassifier(nb, args);
		
	}
	
	
}
