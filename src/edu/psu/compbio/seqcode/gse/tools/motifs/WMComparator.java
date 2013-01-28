package edu.psu.compbio.seqcode.gse.tools.motifs;

import edu.psu.compbio.seqcode.gse.clustering.PairwiseElementMetric;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;

/** compares two weight matrices and returns a score to describe how
 * close they are.  Zero is the minumum score and represents a perfect match.
 *   Worse matches should receive larger values 
 */
public interface WMComparator extends PairwiseElementMetric<WeightMatrix> {

    public double compare(WeightMatrix query, WeightMatrix target);
}
