package edu.psu.compbio.seqcode.projects.akshay.bayesments.features;

import java.io.File;
import java.io.IOException;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.shaun.EventMetaMaker;

public class Sequences {
	
	protected List<Point> locations;
	protected Genome gen;
	protected Config conf;
	protected int winSize;
	
	public Sequences(Config conf) {
		this.conf = conf;
		this.winSize = conf.getSeqWinSize();
		File peaksFile = conf.getPeaksFile();
		try {
			locations = EventMetaMaker.loadPoints(peaksFile);
			this.gen = conf.getGenome();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Sequences(Config conf, List<Point> locations) {
		this.conf =conf;
		this.gen = conf.getGenome();
		this.locations = locations;
		this.winSize = conf.getSeqWinSize();
	}
	
	
	
	
	
}
