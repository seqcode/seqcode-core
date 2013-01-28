/*
 * Created on Mar 10, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import java.util.Iterator;


import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.utils.*;

/**
 * @author tdanford
 * 
 * RegionProber is the interface for objects which turn Regions into lists of
 * Probes.  At the moment, it's an interface which simply extends the Expander
 * interface.  However, it might hold additional, RegionProber-specific methods
 * in the future as well.
 * 
 */
public interface RegionProber<X extends Probe> extends Expander<Region,X> {

    public static class Wrapper<A extends Probe> implements RegionProber<A>, Closeable {
        
        private Expander<Region,A> internal;
        private boolean closed;
        
        public Wrapper(Expander<Region,A> i) { 
            internal = i;
            closed = false;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Expander#execute(java.lang.Object)
         */
        public Iterator<A> execute(Region a) {
            if(closed) { throw new IllegalStateException(); }
            return internal.execute(a);
        } 
        
        public boolean isClosed() { return closed; }
        
        public void close() { 
            if(internal instanceof Closeable) { 
                ((Closeable)internal).close();
            }
            
            closed = true;
        }
    }
}
