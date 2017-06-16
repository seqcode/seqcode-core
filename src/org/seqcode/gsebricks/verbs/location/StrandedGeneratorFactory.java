package org.seqcode.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gsebricks.RegionExpanderFactory;
import org.seqcode.gsebricks.verbs.Expander;

public class StrandedGeneratorFactory implements RegionExpanderFactory<StrandedRegion> {
	String type;

	public StrandedGeneratorFactory() {
		type = "StrandedRegion";
	}

	public StrandedGeneratorFactory(String t) {
		type = t;
	}

	public void setType(String t) {
		type = t;
	}

	public String getType() {
		return type;
	}

	public String getProduct() {
		return "StrandedRegion";
	}

	public Expander<Region, StrandedRegion> getExpander(Genome g) {
		return getExpander(g, type);
	}

	public Expander<Region, StrandedRegion> getExpander(Genome g, String type) {
		if (type == null) {
			throw new NullPointerException("StrandedGenerator must have a type");
		} else {
			return new StrandedRegionGenerator(g, type);
		}
	}
}
