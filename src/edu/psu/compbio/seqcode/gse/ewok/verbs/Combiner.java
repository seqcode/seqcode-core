package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.Iterator;

public interface Combiner<A,B,C> {

    public C execute(A a, B b);

}
