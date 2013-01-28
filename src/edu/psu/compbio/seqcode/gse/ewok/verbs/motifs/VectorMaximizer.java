/*
 * Created on Apr 27, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.motifs;

import java.util.*;

import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.utils.*;

/**
 * @author tdanford
 */
public class VectorMaximizer implements Mapper<Vector<Double>,Pair<Integer,Double>> {

    public VectorMaximizer() {
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper#execute(java.lang.Object)
     */
    public Pair<Integer, Double> execute(Vector<Double> a) {
        int index = -1; 
        double maxValue = 0.0;
        
        for(int i = 0; i < a.size(); i++) { 
            double v = a.get(i);
            if(index == -1 || v > maxValue) { 
                maxValue = v;
                index = i;
            }
        }
        
        return new Pair<Integer,Double>(index, maxValue);
    }

}
