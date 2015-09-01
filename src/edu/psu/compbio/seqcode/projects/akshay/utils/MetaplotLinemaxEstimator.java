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
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
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
	protected int readExt=200;
	
	
	public MetaplotLinemaxEstimator(ExptConfig ef, GenomeConfig gf) {
		econf=ef;
		gconf=gf;
	}
	
	
	//Settors
	public void setPeaks(List<Point> ps){locations = ps;}
	public void setWin(int w){win=w;}
	public void setNBins(int nb){nbins=nb;}
	public void setReadExt(int ext){readExt=ext;}
	
	
	public void execute(){
		BinningParameters params = new BinningParameters(win, nbins);
		ExperimentManager manager = new ExperimentManager(econf);
		List<Sample> sam = manager.getSamples();
		
		List<Double> perBinCounts = new ArrayList<Double>();
		
		for(Point p : locations){
			int left = win/2;
			int right = win-left-1;
			int start  = Math.max(1, p.getLocation()-left);
			int end = Math.min(p.getLocation()+right, p.getGenome().getChromLength(p.getChrom()));
			
			Region query = new Region(p.getGenome(),p.getChrom(),start,end);
			Region extQuery = new Region(p.getGenome(),p.getChrom(),start-readExt >0?start-readExt : 1, end+readExt < p.getGenome().getChromLength(p.getChrom()) ? end+readExt : p.getGenome().getChromLength(p.getChrom()));
			
			double[] array  = new double[params.getNumBins()];
			for(Sample s : sam){
				List<StrandedBaseCount> sbcs = s.getBases(extQuery);
				for(StrandedBaseCount sbc : sbcs){
					SeqHit hit = new SeqHit(p.getGenome(),p.getChrom(),sbc);
					if(readExt>0)
						hit = hit.extendHit(readExt);
					if(hit.overlaps(query)){
						int startOffset = Math.max(0, hit.getStart()-start);
						int endOffset = Math.max(0, Math.min(end, hit.getEnd()-start));
						int tmpEnd = win-startOffset;
						int tmpStart = win-endOffset;
						startOffset = tmpStart;
						endOffset = tmpEnd;
						
						int startbin = params.findBin(startOffset);
						int endbin = params.findBin(endOffset);
						
						addToArray(startbin,endbin,array,1.0);
					}
				}
			}
			for(int i=0; i<array.length; i++){
				perBinCounts.add(array[i]);
			}
		}
		
		Collections.sort(perBinCounts);
		int seventyFiveInd = (int)(perBinCounts.size()*0.75);
		int ninghtyFiveInd = (int)(perBinCounts.size()*0.95);
		int ninghtyNineInd = (int)(perBinCounts.size()*0.99);
		
		System.out.println("75 Percentile PerBinHits:-");
		System.out.println(perBinCounts.get(seventyFiveInd));
		System.out.println("95 Percentile PerBinHits:-");
		System.out.println(perBinCounts.get(ninghtyFiveInd));
		System.out.println("99 Percentile PerBinHits:-");
		System.out.println(perBinCounts.get(ninghtyNineInd));
		
		manager.close();
	}
		
		
		
		
		
		
	
		
	private void addToArray(int i, int j, double[] array, double value) { 
		for(int k = i; k <= j; k++) { 
			array[k] += value;
		}
	}
		

	
	
	public static void main(String[] args){
		GenomeConfig gc = new GenomeConfig(args);
		ExptConfig ec = new ExptConfig(gc.getGenome(),args);
		ArgParser ap = new ArgParser(args);
		
		int w = Args.parseInteger(args, "win", 1000);
		int nb = Args.parseInteger(args, "bins", 100);
		int ext = Args.parseInteger(args, "ext", 200);
				
		List<Point> locs = RegionFileUtilities.loadPeaksFromPeakFile(gc.getGenome(), ap.getKeyValue("peaks"), w);
		
		MetaplotLinemaxEstimator runner = new MetaplotLinemaxEstimator(ec,gc);
		runner.setPeaks(locs);
		runner.setWin(w);
		runner.setNBins(nb);
		runner.setReadExt(ext);
		runner.execute();
		
		
	}
	
	
	

}
