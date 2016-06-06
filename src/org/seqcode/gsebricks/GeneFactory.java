/*
 * Created on Sep 28, 2006
 */
package org.seqcode.gsebricks;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.Expander;

/**
 * @author tdanford
 *
 * <code>GeneFactory</code> returns an Expander than maps Regions to Genes.  The purpose
 * of the Factory is to return an appropriate expander for a given Genome.
 */
public interface GeneFactory {
    public Expander<Region,Gene> getExpander(Genome g);
}
