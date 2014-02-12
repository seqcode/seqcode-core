package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.Sequences;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class EMrunner {
	
	protected EMtrain model;
	protected boolean onlyChrom;
	protected int num_chrom_itrs;
	protected int num_seq_itrs;
	protected int total_itrs;
	protected GenomicLocations chromdata;
	protected Sequences seqdata;
	protected boolean finishedTraining;
	
	protected Config config;
	
	public EMrunner(Config conf, GenomicLocations chromdata, Sequences seqdata) {
		config = conf;
		this.chromdata = chromdata;
		this.seqdata = seqdata;
		this.onlyChrom = conf.runOnlyChrom();
		this.num_chrom_itrs = conf.getNumChromIters();
		this.num_seq_itrs = conf.getNumSeqIters();
		this.total_itrs = conf.getNumItrs();
	}
	
	
	
}
