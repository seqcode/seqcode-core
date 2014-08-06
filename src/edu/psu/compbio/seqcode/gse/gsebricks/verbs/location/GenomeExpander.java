package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import java.util.*;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.CastingMapper;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.ExpanderIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Mapper;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.MapperIterator;

/**
 * This encapsulates the little pattern that we always have to write when we want
 * to get all of some X in an *entire* genome.
 * 
 * We usually take a Genome, throw it through ChromRegionIterator, (sometimes) 
 * map those NamedRegions down to Regions, and then concatenate that Iterator
 * with the Expander of interest.  This just does all that, in its execute method --
 * it is an Expander which takes a Genome, and returns "all of the X's" in that 
 * entire genome.  
 * 
 * @author tdanford
 */
public class GenomeExpander<X> implements Expander<Genome,X> {
	
	private Expander<Region,X> expander;
	private Mapper<NamedRegion,Region> caster;
	
	public GenomeExpander(Expander<Region,X> exp) { 
		expander = exp;
		caster = new CastingMapper<NamedRegion,Region>();
	}

	public Iterator<X> execute(Genome a) {
        ChromRegionIterator chroms = new ChromRegionIterator(a);
        Iterator<Region> rchroms = new MapperIterator<NamedRegion,Region>(caster, chroms);
		return new ExpanderIterator<Region,X>(expander, rchroms);
	}

}
