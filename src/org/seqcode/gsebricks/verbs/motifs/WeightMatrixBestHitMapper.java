package org.seqcode.gsebricks.verbs.motifs;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.Mapper;

public class WeightMatrixBestHitMapper implements Mapper<Region, WeightMatrixHit> {

	private WeightMatrixScorer scorer;

	public WeightMatrixBestHitMapper(WeightMatrixScorer s) {
		scorer = s;
	}

	public WeightMatrixHit execute(Region a) {

		WeightMatrixScoreProfile prof = scorer.execute(a);
		int bestHit = prof.getMaxIndex();
		char bestStrand = prof.getMaxStrand(bestHit);
		double score = prof.getMaxScore(bestHit);

		WeightMatrix m = prof.getMatrix();

		WeightMatrixHit hit = new WeightMatrixHit(a.getGenome(), a.getChrom(), bestHit, bestHit + m.matrix.length - 1,
				score, bestStrand, m);
		return hit;
	}

}
