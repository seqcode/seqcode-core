/*
 * Created on Mar 10, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;

/**
 * @author tdanford
 */
public interface ProbeAnalyzer extends Filter<Probe,Double> {
    
    public static class Index implements ProbeAnalyzer {
        
        private String key;
        private int index;
        
        public Index(String k, int i) { 
            key = k; 
            index = i;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Double execute(Probe a) {
            if(!a.containsKey(key) || index > a.getValue(key).length) { 
                return null;
            }
            return a.getValue(key)[index]; 
        } 
        
    }
    
    public static class BayesPosterior extends Index {
        public BayesPosterior(String k) { 
            super(k, 1);
        }
    }
    
    public static class BayesStrength extends Index {
        public BayesStrength(String k) { 
            super(k, 0);
        }
    }

}
