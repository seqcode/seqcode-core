package org.seqcode.genome.sequence.seqfunctions;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SeqFunctionMask implements SeqFunction {

	SeqFunction function;
	List<Integer> posList; // relative coordinates w.r.t scoring array centers
	Map<Integer, Integer> rel2compactindex; // mapping of relative coord to
											// index in compact representation

	public SeqFunctionMask(SeqFunction function, List<Integer> posList) {
		this.function = function;
		this.posList = posList;

		rel2compactindex = new HashMap<Integer, Integer>();
		int index = 0;
		for (Integer i : posList) {
			rel2compactindex.put(i, index);
			index++;
		}
	}

	// Returns compact/compressed masked representation
	public double[][] score(String seq) throws SeqFunctionException {
		double[][] fScores = function.score(seq);
		double[][] maskScores = new double[function.scoreDimension()][posList.size()];

		for (int d = 0; d < function.scoreDimension(); d++) {
			for (int i = 0; i < seq.length(); i++) {
				int rel = i - (seq.length() / 2);
				if (posList.contains(rel))
					maskScores[d][rel2compactindex.get(rel)] = fScores[d][i];
			}
		}
		return maskScores;
	}

	public List<Integer> getRelPosList() {
		return posList;
	}

	// score dimensionality (per base)
	public int scoreDimension() {
		return function.scoreDimension();
	}

	// offset from initial base position to first score.
	public int scoringOffset() {
		return 0;
	}

	// Window size is the window of sequence that each score is based on
	public int scoreWindowSize() {
		return function.scoreWindowSize();
	}

	// the score is defined between nucleotide positions (i.e. first score
	// defined at offset+0.5bp)
	public boolean isBetweenNucleotides() {
		return function.isBetweenNucleotides();
	}

	// Score max & min
	public double getMaxScore() {
		return function.getMaxScore();
	}

	public double getMinScore() {
		return function.getMinScore();
	}

	// Labels on each dimension
	public String[] dimensionLabels() {
		return function.dimensionLabels();
	}

	// Description of the score
	public String scoreDescription() {
		return function.scoreDescription();
	}
}
