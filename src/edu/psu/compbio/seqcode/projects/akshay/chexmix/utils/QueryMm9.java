package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

public class QueryMm9 extends QueryGenome{

	public QueryMm9(String chr, int midpoint, int range) {
		super(chr, midpoint, range);
	}

	@Override
	public void fillGenomePath() {
		this.genomepath = "/gpfs/home/auk262/group/genomes/mm9";
		
	}

	@Override
	public void fillGenomePath(String path) {
		this.genomepath = path;
	}

}
