package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;

public interface MultiSink<X> { 
	public void init();
	public void consume(String n, X val);
    public void consume(String n, Iterator<X> itr);
	public void finish();
}
