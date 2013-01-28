/*
 * Created on Mar 20, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

/**
 * @author tdanford
 */
public class CastingMapper<A,B> implements Mapper<A,B> { 

    public CastingMapper() {
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(null)
     */
    public B execute(A a) {
        return (B)a;
    }

}
