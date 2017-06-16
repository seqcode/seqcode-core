package org.seqcode.genome.sequence.seqfunctions;

/**
 * Converts the sequence into a 2xL matrix of Strong/Weak base frequencies
 * 
 * @author mahony
 *
 */
public class SWBaseFrequencyFunction implements SeqFunction {

	// Variables
	final int scoreDimension = 2;
	final int scoringOffset = 0;
	final int scoreWindowSize = 1;
	final boolean isBetweenNucs = false;
	final String[] labels = { "S", "W" };
	final String description = "Strong/Weak base frequencies";

	public double[][] score(String seq) throws SeqFunctionException {
		if (seq.length() < scoreWindowSize)
			throw new SeqFunctionException("Sequence too short for SWBaseFrequencyFunction");

		double[][] scores = new double[scoreDimension][seq.length()];
		String seqU = seq.toUpperCase();
		for (int i = 0; i < seqU.length(); i++) {
			char b = seqU.charAt(i);
			for (int x = 0; x < scoreDimension; x++) {
				scores[x][i] = 0;
			}

			switch (b) {
			case 'A':
				scores[1][i] = 1;
				break;
			case 'C':
				scores[0][i] = 1;
				break;
			case 'G':
				scores[0][i] = 1;
				break;
			case 'T':
				scores[1][i] = 1;
				break;
			case 'N':
			case 'Y':
			case 'M':
			case 'R':
			case 'K':
			default:
				scores[0][i] = 0.5;
				scores[1][i] = 0.5;
				break;
			case 'S':
				scores[0][i] = 1;
				break;
			case 'W':
				scores[1][i] = 1;
				break;
			}
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
