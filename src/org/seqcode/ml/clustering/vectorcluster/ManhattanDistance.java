package org.seqcode.ml.clustering.vectorcluster;

import org.seqcode.ml.clustering.PairwiseElementMetric;

public class ManhattanDistance<X extends VectorClusterElement> implements PairwiseElementMetric<X> {
    public ManhattanDistance() {}
    public double evaluate(X e1, X e2) {
        if(e1.dimension() != e2.dimension()) { throw new IllegalArgumentException(); }
        double value = 0.0;
        
        for(int i = 0; i < e1.dimension(); i++) { 
            double s = 0.0;
            if(!e1.isMissingValue(i) && !e2.isMissingValue(i)) { 
                s = (e1.getValue(i) - e2.getValue(i));
            }
            value += Math.abs(s);
        }
			
        return value;
    } 
		
}
