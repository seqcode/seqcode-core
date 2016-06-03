package org.seqcode.projects.akshay.chexmix.utils;

import java.io.*;
import java.util.*;

import org.seqcode.projects.akshay.chexmix.datasets.Seq;

public class QueryHg19 extends QueryGenome{

	public QueryHg19() {
		super();
	}

	@Override
	public void fillGenomePath() {
		this.genomepath = "/gpfs/home/auk262/group/genomes/hg19";
	}

	@Override
	public void fillGenomePath(String path) {
		this.genomepath = path;
	}
}
