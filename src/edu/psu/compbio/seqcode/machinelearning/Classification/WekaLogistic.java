package edu.psu.compbio.seqcode.machinelearning.Classification;

import weka.classifiers.*;
import weka.classifiers.functions.Logistic;



public class WekaLogistic {
		
	public static void main(String[] args){
		Logistic lg = new Logistic();
		AbstractClassifier.runClassifier(lg, args);
	}
	

}
