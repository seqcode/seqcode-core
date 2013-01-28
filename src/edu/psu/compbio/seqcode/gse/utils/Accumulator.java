package edu.psu.compbio.seqcode.gse.utils;

public interface Accumulator<T> {
	public void accumulate(T value);
}
