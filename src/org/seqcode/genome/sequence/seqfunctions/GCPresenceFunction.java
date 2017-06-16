package org.seqcode.genome.sequence.seqfunctions;

/**
 * Indicates presence/absense of *any* G/C bases within a sliding window of
 * length W along a sequence
 * 
 * @author mahony
 *
 */
public class GCPresenceFunction implements SeqFunction {

	// Variables
	final int scoreDimension = 1;
	int scoringOffset = 0;
	int scoreWindowSize = 5;
	boolean isBetweenNucs = false;
	final String[] labels = { "GCpres" };
	final String description = "GC presence";

	public GCPresenceFunction(int W) {
		scoreWindowSize = W;
		if (scoreWindowSize % 2 == 0) {
			scoringOffset = (scoreWindowSize / 2) - 1;
			isBetweenNucs = true;
		} else {
			scoringOffset = scoreWindowSize / 2;
			isBetweenNucs = false;
		}
		String l = String.format("%s_%dbp", "GCpres", scoreWindowSize);
		labels[0] = l;
	}

	public double[][] score(String seq) throws SeqFunctionException {
		if (seq.length() < scoreWindowSize)
			throw new SeqFunctionException("Sequence too short for GCPresenceFunction");

		double[][] scores = new double[scoreDimension][seq.length()];
		double[] composition = new double[seq.length()];
		String seqU = seq.toUpperCase();
		for (int i = 0; i < seqU.length(); i++) {
			char b = seqU.charAt(i);
			scores[0][i] = 0;

			switch (b) {
			case 'A':
				composition[i] = 0;
				break;
			case 'C':
				composition[i] = 1;
				break;
			case 'G':
				composition[i] = 1;
				break;
			case 'T':
				composition[i] = 0;
				break;
			case 'N':
			default:
				composition[i] = 0.5;
				break;
			case 'R':
				composition[i] = 0.5;
				break;
			case 'Y':
				composition[i] = 0.5;
				break;
			case 'M':
				composition[i] = 0.5;
				break;
			case 'K':
				composition[i] = 0.5;
				break;
			case 'S':
				composition[i] = 1;
				break;
			case 'W':
				composition[i] = 0;
				break;
			}
		}

		for (int i = 0; i < composition.length - scoreWindowSize + 1; i++) {
			double score = 0;
			for (int w = 0; w < scoreWindowSize; w++) {
				score += composition[i + w];
			}
			if (score > 0)
				score = 1;
			scores[0][i + scoringOffset] = score;
		}

		return scores;
	}

	public int scoreDimension() {
		return scoreDimension;
	}

	public int scoringOffset() {
		return scoringOffset;
	}

	public int scoreWindowSize() {
		return scoreWindowSize;
	}

	public boolean isBetweenNucleotides() {
		return isBetweenNucs;
	}

	public String[] dimensionLabels() {
		return labels;
	}

	public String scoreDescription() {
		return description;
	}

	public double getMaxScore() {
		return 1.0;
	}

	public double getMinScore() {
		return 0.0;
	}
}
