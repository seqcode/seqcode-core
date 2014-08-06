package edu.psu.compbio.seqcode.gse.gsebricks.types;

import java.util.Collection;

public interface SelfDescribingOutput<X> extends SelfDescribingVerb { 
	public Collection<X> getValues();
}
