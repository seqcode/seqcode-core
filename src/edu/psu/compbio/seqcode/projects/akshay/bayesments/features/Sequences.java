package edu.psu.compbio.seqcode.projects.akshay.bayesments.features;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.shaun.EventMetaMaker;

public class Sequences {
	
	protected List<Point> locations;
	protected Genome gen;
	protected Config conf;
	protected int winSize;
	protected String[] sequences;
	protected double Xs[][];
	protected List<WeightMatrix> motifs_log_odds;
	
	public Sequences(Config conf, List<WeightMatrix> wm) {
		this.conf = conf;
		this.winSize = conf.getSeqWinSize();
		this.motifs_log_odds = wm;
		File peaksFile = conf.getPeaksFile();
		try {
			locations = EventMetaMaker.loadPoints(peaksFile);
			this.gen = conf.getGenome();
			this.fetchSequences();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Sequences(Config conf, List<Point> locations, List<WeightMatrix> wm) {
		this.conf =conf;
		this.gen = conf.getGenome();
		this.locations = locations;
		this.winSize = conf.getSeqWinSize();
		this.fetchSequences();
		this.motifs_log_odds = wm;
	}
	
	private void fetchSequences(){
		Region[] regions = new Region[locations.size()];
		int countPoint = 0;
		for(Point p : locations){
			regions[countPoint] =  new Region(this.gen,p.getChrom()
					,p.getLocation()-winSize,p.getLocation()+winSize);
		}
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(true);
		seqgen.useLocalFiles(true);
		seqgen.setGenomePath(conf.getGenomeSeqPath());
		sequences = new String[locations.size()];
		for(int i=0; i< regions.length; i++){
			sequences[i] =seqgen.execute(regions[i]);
			sequences[i].toLowerCase();
		}
	}
	
	public void setXc(){
		Xs =  new double[locations.size()][motifs_log_odds.size()];
		
	}
	
	private double calculateLogOddinRegion(WeightMatrix wm, String seq){
		int motif_width = wm.length();
		int seq_length = seq.length();
		double total_score=0.0;
		for(int i=0; i< seq_length-motif_width; i++){
			String scan_seq = seq.substring(i, i+motif_width);
			double score =1.0;
			for(int j=0; j<motif_width; j++){
				String letter = scan_seq.substring(j,j+1);
				
			}
		}
		return total_score;
	}
	
	
	
	
	
}
