package org.seqcode.data.seqdata;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;

public class SeqAnalysisResult extends Region {
	public Integer position;
	public Double foregroundReadCount, backgroundReadCount, strength, shape, pvalue, foldEnrichment;

	public SeqAnalysisResult(Genome g, String chrom, int start, int end, Integer position, Double fgcount,
			Double bgcount, Double strength, Double shape, Double pvalue, Double foldEnrichment) {
		super(g, chrom, start, end);
		this.position = position;
		this.foregroundReadCount = fgcount;
		this.backgroundReadCount = bgcount;
		this.strength = strength;
		this.shape = shape;
		this.pvalue = pvalue;
		this.foldEnrichment = foldEnrichment;
	}

	public Integer getPosition() {
		return position;
	}

	public Double getFG() {
		return foregroundReadCount;
	}

	public Double getBG() {
		return backgroundReadCount;
	}

	public Double getStrength() {
		return strength;
	}

	public Double getShape() {
		return shape;
	}

	public Double getPValue() {
		return pvalue;
	}

	public Double getFoldEnrichment() {
		return foldEnrichment;
	}
}