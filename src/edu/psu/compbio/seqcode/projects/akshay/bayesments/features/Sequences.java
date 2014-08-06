package edu.psu.compbio.seqcode.projects.akshay.bayesments.features;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.EventMetaMaker;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.CoutPlotter;

public class Sequences {
	
	protected List<Point> locations;
	protected Genome gen;
	protected Config conf;
	protected int winSize;
	protected String[] sequences;
	protected double Xs[][];
	protected List<WeightMatrix> motifs_log_odds;
	
	public Sequences(Config conf) {
		this.conf = conf;
		this.winSize = conf.getSeqWinSize();
		File peaksFile = conf.getPeaksFile();
		try {
			this.gen = conf.getGenome();
			locations = EventMetaMaker.loadPoints(peaksFile, this.gen);
			this.fetchSequences();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Sequences(Config conf, List<Point> locations) {
		this.conf =conf;
		this.gen = conf.getGenome();
		this.locations = locations;
		this.winSize = conf.getSeqWinSize();
		this.fetchSequences();
	}
	
	private void fetchSequences() {
		Region[] regions = new Region[locations.size()];
		int countPoint = 0;
		for(Point p : locations){
			regions[countPoint] = p.expand(winSize);
			countPoint++;
			
		}
		
		
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useLocalFiles(true);
		seqgen.setGenomePath(conf.getGenomeSeqPath());
		seqgen.useCache(true);
		sequences = new String[locations.size()];
		for(int i=0; i< regions.length; i++){
			sequences[i] =seqgen.execute(regions[i]);
			sequences[i].toLowerCase();
			
		}
	}
	
	//setters
	public void setXc(){
		Xs =  new double[locations.size()][motifs_log_odds.size()];
		for(int i=0; i<locations.size(); i++){
			for(int j=0; j<motifs_log_odds.size(); j++){
				Xs[i][j] = this.calculateLogOddinRegion(this.motifs_log_odds.get(j), this.sequences[i]);
			}
		}
	}
	
	public void setMotifs(List<WeightMatrix> wm){this.motifs_log_odds = wm;}
	
	public void plotSeqScores(){
		for(int m=0; m<	this.motifs_log_odds.size()	; m++){
			float[] counts = new float[locations.size()];
			for(int i=0; i<locations.size(); i++){
				counts[i] = (float)Xs[i][m];
			}
			String name_tag = "Motif_"+Integer.toString(m);
			CoutPlotter cp = new CoutPlotter(counts, conf, name_tag);
			cp.plot();
		}
	}
	
	private double calculateLogOddinRegion(WeightMatrix wm, String seq){
		int motif_width = wm.length();
		int seq_length = seq.length();
		//given string
		double total_score=-10000.0;
		for(int i=0; i< seq_length-motif_width+1; i++){
			double score =0.0;
			for(int j=0; j<motif_width; j++){
				score = score+wm.matrix[j][seq.charAt(i+j)];
			}
			total_score = (total_score < score ? score : total_score);
		}
		//rev complement 
		String rev_seq = SequenceUtils.reverseComplement(seq);
		for(int i=0; i< seq_length-motif_width+1; i++){
			double score =0.0;
			for(int j=0; j<motif_width; j++){
				score = score+wm.matrix[j][rev_seq.charAt(i+j)];
			}
			total_score = (total_score < score ? score : total_score);
		}
		if(conf.doSeqAshin()){
			total_score = Math.log(total_score+ Math.sqrt(Math.pow(total_score, 2.0)+1.0));
		}
		return total_score;
	}
	
	//Accessors
	
	public double[][] getXs(){return this.Xs;}
	public int getNumMotifs(){return this.motifs_log_odds.size();}
	// i index starts from 0
	public String getIthSeq(int i){return this.sequences[i];}
	
	public static void main(String[] args){
		double s = -2.0;
		double tt=10;
		
		
		System.out.println(Math.log(tt));
	}
	
	
}
