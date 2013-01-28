package edu.psu.compbio.seqcode.gse.ewok.verbs;

public interface BiCombiner<X,Y,Z> {
	public Z execute(X a, Y b);
}
