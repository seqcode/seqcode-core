/*
 * Created on Mar 10, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;

/**
 * @author tdanford
 */
public interface PeakFinder<X extends Probe> extends Distiller<X,BindingEvent> {
    
    public static class Wrapper<A extends Probe> implements PeakFinder<A> {
        
        private Distiller<A,BindingEvent> dist;
        
        public Wrapper(Distiller<A,BindingEvent> d) { 
            dist = d;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#execute(java.lang.Object)
         */
        public BindingEvent execute(A a) {
            return dist.execute(a);
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#getCurrent()
         */
        public BindingEvent getCurrent() {
            return dist.getCurrent();
        }

        public void reset() {
            dist.reset();
        } 
        
    }
}
