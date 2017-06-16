package org.seqcode.ml.classification;

import weka.classifiers.*;
import weka.classifiers.functions.Logistic;

public class WekaLogistic {

	public static void main(String[] args) {
		Logistic lg = new Logistic();
		AbstractClassifier.runClassifier(lg, args);
	}

}
