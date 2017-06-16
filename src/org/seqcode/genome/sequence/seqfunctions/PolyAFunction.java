package org.seqcode.genome.sequence.seqfunctions;

/**
 * Scores poly-A (or poly-T) instances of length W along a sequence
 * 
 * @author mahony
 *
 */
public class PolyAFunction implements SeqFunction {

	// Variables
	final int scoreDimension = 1;
	int scoringOffset = 0;
	int scoreWindowSize = 5;
	boolean isBetweenNucs = false;
	final String[] labels = { "PolyA" };
	final String description = "PolyA";

	public PolyAFunction(int W) {
		scoreWindowSize = W;
		if (scoreWindowSize % 2 == 0) {
			scoringOffset = (scoreWindowSize / 2) - 1;
			isBetweenNucs = true;
		} else {
			scoringOffset = scoreWindowSize / 2;
			isBetweenNucs = false;
		}
		String l = String.format("%s_%dbp", "PolyA", scoreWindowSize);
		labels[0] = l;
	}

	public double[][] score(String seq) throws SeqFunctionException {
		if (seq.length() < scoreWindowSize)
			throw new SeqFunctionException("Sequence too short for PolyAFunction");

		double[][] scores = new double[scoreDimension][seq.length()];
		double[] polyA = new double[seq.length()];
		double[] polyT = new double[seq.length()];
		String seqU = seq.toUpperCase();
		for (int i = 0; i < seqU.length(); i++) {
			char b = seqU.charAt(i);
			scores[0][i] = 0;

			switch (b) {
			case 'A':
				polyA[i] = 1;
				break;
			case 'T':
				polyT[i] = 1;
				break;
			default:
				polyA[i] = 0;
				polyT[i] = 0;
				break;
			}
		}

		for (int i = 0; i < polyA.length - scoreWindowSize + 1; i++) {
			int A = 0, T = 0;
			for (int w = 0; w < scoreWindowSize; w++) {
				A += polyA[i + w];
				T += polyT[i + w];
			}
			if (A == scoreWindowSize || T == scoreWindowSize)
				scores[0][i + scoringOffset] = 1;
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
