/*
 * Created on Sep 29, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.PeakCallerFactory;
import edu.psu.compbio.seqcode.gse.ewok.verbs.probers.ChipChipImmediateProbeGenerator;

/**
 * @author tdanford
 */
public class ThresholdPeakCallerFactory implements PeakCallerFactory {
    
    public ThresholdPeakCallerFactory() {
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.PeakCallerFactory#createCaller(edu.psu.compbio.seqcode.gse.Bio.Genome, edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator)
     */
    public PeakCaller createCaller(Genome g, ExptLocator loc, Object args) {
        ChipChipLocator aloc = (ChipChipLocator)loc;
        ChipChipImmediateProbeGenerator gen = new ChipChipImmediateProbeGenerator(g, aloc);
        RegionProber<Probe> prober = new RegionProber.Wrapper<Probe>(gen);
        
        String ename = aloc.name + "," + aloc.version;
        double thresh = args != null ? (Double)args : 2.0;
        ThresholdPeakFinder<Probe> finder = new ThresholdPeakFinder<Probe>(ename, thresh);
        
        PeakCaller caller = new PeakCaller.FromFinder<Probe>(prober, finder);
        return caller;
    }

}
