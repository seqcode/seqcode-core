package edu.psu.compbio.seqcode.gse.gsebricks.verbs;

public interface BiCombiner<X,Y,Z> {
	public Z execute(X a, Y b);
}
