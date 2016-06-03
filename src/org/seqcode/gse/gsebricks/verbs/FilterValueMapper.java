/*
 * Created on Mar 10, 2006
 */
package org.seqcode.gse.gsebricks.verbs;

import org.seqcode.gse.utils.*;

/**
 * @author tdanford
 */
public class FilterValueMapper<X,Y> implements Mapper<X,Pair<X,Integer>> {
    
    private Filter<X,Y> filter;

    public FilterValueMapper(Filter<X,Y> f) {
        filter = f;
    }

    /* (non-Javadoc)
     * @see org.seqcode.gse.gsebricks.verbs.Mapper#execute(java.lang.Object)
     */
    public Pair<X, Integer> execute(X a) {
        Y ret = filter.execute(a);
        if(ret != null) { 
            return new Pair<X,Integer>(a, 1); 
        } else { 
            return new Pair<X,Integer>(a, 0);
        }
    }

}
