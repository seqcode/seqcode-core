package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.projects.seed.stats.FeatureStatistics;


/**
 * Simple class to count the number of tags mapping to a given set of regions for a given list of sequencing experiments
 * @author akshaykakumanu
 *
 */
public class SequencingExptRegionsCounter {
	
	ExptConfig econf;
	GenomeConfig gconf;
	
	double minSigCtrlFoldDifference;
	
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
	
	
	public void printInputEnrichment(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		double[] enrichment = new double[regions.size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		
		for(ExperimentCondition ec: manager.getConditions()){
			header.append(ec.getName());
			header.append("\t");
			double[] pooledSigCounts= new double[regions.size()];
			double[] pooledCtrlCounts = new double[regions.size()];
			
			for(Sample sigs : ec.getSignalSamples()){
				int RegionCount=0;
				for(Region r : regions){
					pooledSigCounts[RegionCount] = pooledSigCounts[RegionCount] + sigs.countHits(r);
					RegionCount++;
				}
			}
			
			if(ec.getControlSamples().size() != 0 ){
				for(Sample ctrls : ec.getControlSamples()){
					int RegionCount=0;
					for(Region r : regions){
						pooledCtrlCounts[RegionCount] = pooledCtrlCounts[RegionCount] + ctrls.countHits(r);
						RegionCount++;
					}
				}
			}else{
				int RegionCount=0;
				for(Region r : regions){
					pooledCtrlCounts[RegionCount] = ec.getTotalSignalCount()*r.getWidth()/econf.getMappableGenomeLength();
					RegionCount++;
				}
			}
			// Scale the control counts
			for(Region r : regions){
				int RegionCount=0;
				pooledCtrlCounts[RegionCount] = pooledCtrlCounts[RegionCount]*ec.getPooledSampleControlScaling();
				RegionCount++;
			}
			
			FeatureStatistics stats = new FeatureStatistics();
			// Do bionomial testing 
			for(Region r : regions){
				int RegionCount=0;
				double pval = stats.binomialPValue(pooledCtrlCounts[RegionCount], pooledCtrlCounts[RegionCount]+
						pooledSigCounts[RegionCount],minSigCtrlFoldDifference);
				enrichment[RegionCount] = -1*Math.log10(pval);
				RegionCount++;
			}
		}
		
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		
		StringBuilder enrichmentSb = new StringBuilder();
		
		for(int r=0; r<regions.size(); r++){
			enrichmentSb.append(regions.get(r).getLocationString());
			enrichmentSb.append("\t");
			enrichmentSb.append(enrichment[r]);
			enrichmentSb.append("\n");
		}
		System.out.println(enrichmentSb.toString());
		
		manager.close();
	}
	
	
	
	public void setRegions(List<Region> regs){regions = regs;}
	public void serMinSigCntrlFoldDiff(double min){minSigCtrlFoldDifference = min;}
	
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		
		ArgParser ap = new ArgParser(args);
		
		SequencingExptRegionsCounter counter = new SequencingExptRegionsCounter(econ, gcon);
		
		int win = Args.parseInteger(args, "win", 200);
		
		double minSigCtrlFoldDifference = Args.parseDouble(args,"minfolddiff",1);
		
		if (!ap.hasKey("refTSSs") && !ap.hasKey("peaks") && !ap.hasKey("regions")){
			System.err.println("Provied either a peaks file or a refTSSs file or a regions !!!");
			System.exit(1);
		}
		
		List<Region> reg = new ArrayList<Region>();
		
		if(ap.hasKey("peaks") || ap.hasKey("refTSSs")){
			@SuppressWarnings("unchecked")
			List<Point> points = (List<Point>) (ap.hasKey("refTSSs") ? RegionFileUtilities.loadStrandedPointFromRefTssFile(gcon.getGenome(), ap.getKeyValue("refTSSs"))
					: RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), ap.getKeyValue("peaks"), win));
			for(Point p: points){
				reg.add(p.expand(win/2));
			}
		}else{
			reg = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), ap.getKeyValue("regions"), -1);
		}
		
		
		counter.setRegions(reg);
		counter.serMinSigCntrlFoldDiff(minSigCtrlFoldDifference);
		if(ap.hasKey("counts")){
			counter.printCounts();
		}else if(ap.hasKey("enrichment")){
			counter.printInputEnrichment();
		}
		
		
		
	}

	
	

}
