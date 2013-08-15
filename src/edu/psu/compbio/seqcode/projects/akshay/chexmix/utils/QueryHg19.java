package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Seq;
import java.io.*;
import java.util.*;

public class QueryHg19 extends QueryGenome{

	public QueryHg19(String chr, int midpoint, int range) {
		super(chr, midpoint, range);
	}

	@Override
	public void fillGenomePath() {
		this.genomepath = "/gpfs/home/auk262/group/genomes/hg19";
	}

	@Override
	public void fillGenomePath(String path) {
		this.genomepath = path;
	}
	
	public static void main(String[] args){
		String con = null;
		System.out.println(con+"akshay");
	}

}
