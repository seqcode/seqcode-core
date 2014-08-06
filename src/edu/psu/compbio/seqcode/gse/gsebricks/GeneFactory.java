/*
 * Created on Sep 28, 2006
 */
package edu.psu.compbio.seqcode.gse.gsebricks;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

/**
 * @author tdanford
 *
 * <code>GeneFactory</code> returns an Expander than maps Regions to Genes.  The purpose
 * of the Factory is to return an appropriate expander for a given Genome.
 */
public interface GeneFactory {
    public Expander<Region,Gene> getExpander(Genome g);
}
