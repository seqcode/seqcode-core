package org.seqcode.gsebricks.verbs;

public interface BiCombiner<X,Y,Z> {
	public Z execute(X a, Y b);
}
