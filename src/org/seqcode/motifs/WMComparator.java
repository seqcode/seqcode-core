package org.seqcode.motifs;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.ml.clustering.PairwiseElementMetric;

/**
 * compares two weight matrices and returns a score to describe how close they
 * are. Zero is the minumum score and represents a perfect match. Worse matches
 * should receive larger values
 */
public interface WMComparator extends PairwiseElementMetric<WeightMatrix> {

	public double compare(WeightMatrix query, WeightMatrix target);
}
