package org.seqcode.gsebricks.verbs.motifs;

import java.util.*;

import org.seqcode.data.motifdb.*;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.Expander;
import org.seqcode.gsebricks.verbs.Filter;
import org.seqcode.gsebricks.verbs.Mapper;

public class WeightMatrixFilter implements Filter<Region, Region> {

	private Vector<WeightMatrix> matrices;
	private Vector<Float> cutoffs;
	private SequenceGenerator seqgen;

	public WeightMatrixFilter(WeightMatrix m, float c) {
		matrices = new Vector<WeightMatrix>();
		cutoffs = new Vector<Float>();
		matrices.add(m);
		cutoffs.add(c);
		seqgen = new SequenceGenerator();
	}

	public WeightMatrixFilter() {
		matrices = new Vector<WeightMatrix>();
		cutoffs = new Vector<Float>();
		seqgen = new SequenceGenerator();
	}

	public void addWeightMatrix(WeightMatrix m, float c) {
		matrices.add(m);
		cutoffs.add(c);
	}

	public Region execute(Region r) {
		String seq = seqgen.execute(r);
		float[] scores = null;
		int s = r.getStart();

		try {
			for (int k = 0; k < matrices.size(); k++) {
				WeightMatrix matrix = matrices.get(k);
				int width = matrix.matrix.length;
				float cutoff = cutoffs.get(k);

				if (score(matrix, seq.toCharArray(), cutoff)) {
					return r;
				}
			}

			seq = SequenceUtils.reverseComplement(seq);

			for (int k = 0; k < matrices.size(); k++) {
				WeightMatrix matrix = matrices.get(k);
				int width = matrix.matrix.length;
				float cutoff = cutoffs.get(k);
				if (score(matrix, seq.toCharArray(), cutoff)) {
					return r;
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			e.printStackTrace(System.err);
		}

		return null;
	}

	public boolean score(WeightMatrix matrix, char[] sequence, float thresh) {
		/* scan through the sequence */
		int length = matrix.length();
		for (int i = 0; i < sequence.length - length; i++) {
			float score = (float) 0.0;
			for (int j = 0; j < length; j++) {
				score += matrix.matrix[j][sequence[i + j]];
			}
			if (score >= thresh) {
				return true;
			}
		}
		return false;
	}

}
