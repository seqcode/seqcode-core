/*
 * Created on Mar 11, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

/**
 * @author tdanford
 */
public class ConditionalFilter<X,Y> implements Filter<X,X> {
    
    private Filter<X,Y> internal;

    public ConditionalFilter(Filter<X,Y> intern) {
        internal = intern;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(null)
     */
    public X execute(X a) {
        Y val = internal.execute(a);
        if(val != null) { return a; }
        return null;
    }

}
