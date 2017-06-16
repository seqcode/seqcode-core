package org.seqcode.gsebricks.verbs;

public interface Combiner<A, B, C> {

	public C execute(A a, B b);

}
