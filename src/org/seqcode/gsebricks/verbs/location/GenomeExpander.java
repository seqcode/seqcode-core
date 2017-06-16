package org.seqcode.gsebricks.verbs.location;

import java.util.*;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.CastingMapper;
import org.seqcode.gsebricks.verbs.Expander;
import org.seqcode.gsebricks.verbs.ExpanderIterator;
import org.seqcode.gsebricks.verbs.Mapper;
import org.seqcode.gsebricks.verbs.MapperIterator;

/**
 * This encapsulates the little pattern that we always have to write when we
 * want to get all of some X in an *entire* genome.
 * 
 * We usually take a Genome, throw it through ChromRegionIterator, (sometimes)
 * map those NamedRegions down to Regions, and then concatenate that Iterator
 * with the Expander of interest. This just does all that, in its execute method
 * -- it is an Expander which takes a Genome, and returns "all of the X's" in
 * that entire genome.
 * 
 * @author tdanford
 */
public class GenomeExpander<X> implements Expander<Genome, X> {

	private Expander<Region, X> expander;
	private Mapper<NamedRegion, Region> caster;

	public GenomeExpander(Expander<Region, X> exp) {
		expander = exp;
		caster = new CastingMapper<NamedRegion, Region>();
	}

	public Iterator<X> execute(Genome a) {
		ChromRegionIterator chroms = new ChromRegionIterator(a);
		Iterator<Region> rchroms = new MapperIterator<NamedRegion, Region>(caster, chroms);
		return new ExpanderIterator<Region, X>(expander, rchroms);
	}

}
