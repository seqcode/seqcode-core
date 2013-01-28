/*
 * Created on Sep 29, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok;

import edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.PeakCaller;

/**
 * @author tdanford
 */
public interface PeakCallerFactory {
    public PeakCaller createCaller(Genome g, ExptLocator loc, Object args);
}
