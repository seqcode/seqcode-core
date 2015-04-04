package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;


/**
 * Simple class to count the number of tags mapping to a given set of regions for a given list of sequencing experiments
 * @author akshaykakumanu
 *
 */
public class SequencingExptRegionsCounter {
	
	ExptConfig econf;
	GenomeConfig gconf;
	
	List<Region> regions;
	
	public SequencingExptRegionsCounter(ExptConfig econfig, GenomeConfig gconfig) {
		econf = econfig;
		gconf = gconfig;
	}
	
	public void printCounts(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		float[][] Counts = new float[regions.size()][manager.getSamples().size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		int SampleCount = 0;
		for(Sample sample : manager.getSamples()){
			header.append(sample.getName());
			header.append("\t");
			int RegionCount = 0;
			for(Region r : regions){
				Counts[RegionCount][SampleCount] = sample.countHits(r);
				RegionCount++;
			}
			SampleCount++;
		}
		
		
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		
		StringBuilder CountsSb = new StringBuilder();
		for(int r=0; r<regions.size(); r++){
			CountsSb.append(regions.get(r).getLocationString());
			CountsSb.append("\t");
			for(int c=0; c<manager.getSamples().size(); c++){
				CountsSb.append(Counts[r][c]);
				CountsSb.append("\t");
			}
			CountsSb.deleteCharAt(CountsSb.length()-1);
			CountsSb.append("\n");
		}
		System.out.println(CountsSb.toString());
		
		manager.close();
	}
	
	
	public void setRegions(List<Region> regs){regions = regs;}
	
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		
		ArgParser ap = new ArgParser(args);
		
		SequencingExptRegionsCounter counter = new SequencingExptRegionsCounter(econ, gcon);
		
		int win = Args.parseInteger(args, "win", 200);
		
		if (!ap.hasKey("refTSSs") && !ap.hasKey("peaks")){
			System.err.println("Provied either a peaks file or a refTSSs file !!!");
			System.exit(1);
		}
		@SuppressWarnings("unchecked")
		List<Point> points = (List<Point>) (ap.hasKey("refTSSs") ? RegionFileUtilities.loadStrandedPointFromRefTssFile(gcon.getGenome(), ap.getKeyValue("refTSSs"))
				: RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), ap.getKeyValue("peaks"), win));
		
		List<Region> reg = new ArrayList<Region>();
		
		for(Point p: points){
			reg.add(p.expand(win/2));
		}
		
		counter.setRegions(reg);
		counter.printCounts();
		
		
	}

	
	

}
