/*
 * Created on Sep 28, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;

/**
 * @author tdanford
 *
 * <code>GeneFactory</code> returns an Expander than maps Regions to Genes.  The purpose
 * of the Factory is to return an appropriate expander for a given Genome.
 */
public interface GeneFactory {
    public Expander<Region,Gene> getExpander(Genome g);
}
