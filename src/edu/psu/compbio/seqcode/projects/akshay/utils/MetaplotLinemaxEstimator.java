package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;


import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.BinningParameters;

public class MetaplotLinemaxEstimator {
	
	protected ExptConfig econf;
	protected GenomeConfig gconf;
	
	protected List<Point> locations;
	protected int win;
	protected int nbins;
	
	
	public MetaplotLinemaxEstimator(ExptConfig ef, GenomeConfig gf) {
		econf=ef;
		gconf=gf;
	}
	
	
	//Settors
	public void setPeaks(List<Point> ps){locations = ps;}
	public void setWin(int w){win=w;}
	public void setNBins(int nb){nbins=nb;}
	
	
	public void execute(){
		BinningParameters params = new BinningParameters(win, nbins);
		List<Integer> perbaseHits = new ArrayList<Integer>();
		List<Region> rs = new ArrayList<Region>();
		for(Point p : locations){
			rs.add(p.expand(win/2));
		}
		ExperimentManager manager = new ExperimentManager(econf);
		List<Sample> sam = manager.getSamples();
		for(Sample s : sam){
			for(Region r : rs){
				for(int b=0; b<params.getNumBins(); b++){
					int start = r.getStart()+b*params.getBinSize();
					int end = r.getStart()+(b+1)*params.getBinSize() < r.getGenome().getChromLength(r.getChrom()) ? r.getStart()+(b+1)*params.getBinSize() : r.getGenome().getChromLength(r.getChrom());  
					Region currBin = new Region(r.getGenome(),r.getChrom(),start,end);
					int counts = (int)s.countHits(currBin);
					perbaseHits.add(counts);
				}
			}
		}
		
		Collections.sort(perbaseHits);
		int seventyFiveInd = (int)(perbaseHits.size()*0.75);
		int ninghtyFiveInd = (int)(perbaseHits.size()*0.995);
		System.out.println("75 Percentile PerBaseHits: -");
		System.out.println(perbaseHits.get(seventyFiveInd));
		System.out.println("95 Percentile PerBaseHits: -");
		System.out.println(perbaseHits.get(ninghtyFiveInd));
		manager.close();
	}
	
	
	public static void main(String[] args){
		GenomeConfig gc = new GenomeConfig(args);
		ExptConfig ec = new ExptConfig(gc.getGenome(),args);
		ArgParser ap = new ArgParser(args);
		
		int w = Args.parseInteger(args, "win", 1000);
		int nb = Args.parseInteger(args, "bins", 100);
				
		List<Point> locs = RegionFileUtilities.loadPeaksFromPeakFile(gc.getGenome(), ap.getKeyValue("peaks"), w);
		
		MetaplotLinemaxEstimator runner = new MetaplotLinemaxEstimator(ec,gc);
		runner.setPeaks(locs);
		runner.setWin(w);
		runner.setNBins(nb);
		runner.execute();
		
		
	}
	
	
	

}
