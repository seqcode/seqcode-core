package edu.psu.compbio.seqcode.gse.gsebricks.verbs;

import java.util.Iterator;

public interface CombinerExpander<A,B,C> {

    public Iterator<C> execute(A a, B b);

}
