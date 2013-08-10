package edu.psu.compbio.seqcode.projects.akshay.chexmix.functions;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Seq;

public class QueryHg19 extends QueryGenome{

	public QueryHg19(String chr, int midpoint) {
		super(chr, midpoint);
	}

	@Override
	public void fillGenomePath() {
		this.genomepath = "/gpfs/home/auk262/group/genomes/hg19";
	}

	@Override
	public Seq getSeq() {
		String currdir = System.getProperty("user.dir");
		
		return null;
	}

	@Override
	public void fillGenomePath(String path) {
		this.genomepath = path;
	}
	
	public static void main(String[] args){
		String currdir =  System.getProperty("user.dir");
		System.out.println(currdir);
	}
	

}
