package edu.psu.compbio.seqcode.gse.clustering.vectorcluster;

import edu.psu.compbio.seqcode.gse.clustering.PairwiseElementMetric;

public class EuclideanDistance<X extends VectorClusterElement> implements PairwiseElementMetric<X> {
    public EuclideanDistance() {}
    public double evaluate(X e1, X e2) {
        if(e1.dimension() != e2.dimension()) { throw new IllegalArgumentException(); }
        double value = 0.0;
        
        for(int i = 0; i < e1.dimension(); i++) { 
            double s = 0.0;
            if(!e1.isMissingValue(i) && !e2.isMissingValue(i)) { 
                s = (e1.getValue(i) - e2.getValue(i));
            }
            value += (s * s);
        }
			
        value = Math.sqrt(value);
        return value;
    } 
		
}
