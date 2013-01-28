package edu.psu.compbio.seqcode.gse.ewok.nouns;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class SimpleDomain extends Region {

	public SimpleDomain(Genome g, String c, int start, int end) {
		super(g, c, start, end);
	}	
	
	public SimpleDomain(Region r) { 
		super(r);
	}
	
	public BindingEvent getBindingEvent() { 
		return new BindingEvent(getGenome(), getChrom(), getStart(), getEnd(), 1.0, 1.0, "SimpleDomain");
	}
}
