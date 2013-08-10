package edu.psu.compbio.seqcode.projects.akshay.chexmix.functions;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;

public abstract class QueryGenome {
	public String genomepath;
	private String chr;
	private int midpoint;
	public abstract void fillGenomePath();
	public abstract void fillGenomePath(String path);
	public QueryGenome(String chr, int midpoint) {
		this.chr = chr;
		this.midpoint = midpoint;
	}
	public abstract Seq getSeq();

}
