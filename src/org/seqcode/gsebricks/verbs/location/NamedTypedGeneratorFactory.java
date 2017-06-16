/*
 * Created on Sep 28, 2006
 */
package org.seqcode.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedTypedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.RegionExpanderFactory;
import org.seqcode.gsebricks.verbs.Expander;

public class NamedTypedGeneratorFactory implements RegionExpanderFactory<NamedTypedRegion> {
	String type;

	public NamedTypedGeneratorFactory() {
		type = "NamedTypedRegion";
	}

	public NamedTypedGeneratorFactory(String t) {
		type = t;
	}

	public void setType(String t) {
		type = t;
	}

	public String getType() {
		return type;
	}

	public String getProduct() {
		return "NamedTypedRegion";
	}

	public Expander<Region, NamedTypedRegion> getExpander(Genome g) {
		return getExpander(g, type);
	}

	public Expander<Region, NamedTypedRegion> getExpander(Genome g, String type) {
		if (type == null) {
			throw new NullPointerException("NamedTypedGenerator must have a type");
		} else {
			return new NamedTypedGenerator(g, type);
		}
	}
}
