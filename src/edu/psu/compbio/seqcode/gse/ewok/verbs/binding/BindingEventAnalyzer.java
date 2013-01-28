/*
 * Created on Mar 10, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;

/**
 * @author tdanford
 */
public interface BindingEventAnalyzer extends Mapper<BindingEvent,Double> {
    
    public static class Size implements BindingEventAnalyzer {
        
        public Size() {}

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Double execute(BindingEvent a) {
            return a.getSize();
        } 
        
    }
    
    public static class Confidence implements BindingEventAnalyzer {
        
        public Confidence() {}

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Double execute(BindingEvent a) {
            return a.getConf();
        } 
        
    }
}
