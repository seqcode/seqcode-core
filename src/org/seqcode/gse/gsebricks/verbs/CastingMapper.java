/*
 * Created on Mar 20, 2006
 */
package org.seqcode.gse.gsebricks.verbs;

/**
 * @author tdanford
 */
public class CastingMapper<A,B> implements Mapper<A,B> { 

    public CastingMapper() {
    }

    /* (non-Javadoc)
     * @see org.seqcode.gse.gsebricks.verbs.Filter#execute(null)
     */
    public B execute(A a) {
        return (B)a;
    }

}
