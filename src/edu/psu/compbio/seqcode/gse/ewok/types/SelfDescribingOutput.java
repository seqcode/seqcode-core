package edu.psu.compbio.seqcode.gse.ewok.types;

import java.util.Collection;

public interface SelfDescribingOutput<X> extends SelfDescribingVerb { 
	public Collection<X> getValues();
}
