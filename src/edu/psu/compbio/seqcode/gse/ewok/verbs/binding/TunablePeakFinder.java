/*
 * Created on Mar 10, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;

/**
 * @author tdanford
 */
public interface TunablePeakFinder<X extends Probe> extends PeakFinder<X> {
    public double getMaxParameter();
    public double getMinParameter();
    public double getCurrentParameter();
    public void setParameter(double p);
}
