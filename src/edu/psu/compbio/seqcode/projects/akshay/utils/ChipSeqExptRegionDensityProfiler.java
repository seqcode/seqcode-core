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
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class ChipSeqExptRegionDensityProfiler {
	
	private GenomeConfig gconf;
	private ExptConfig econf;
	List<Region> regions = new ArrayList<Region>();
	
	public ChipSeqExptRegionDensityProfiler(ExptConfig econ, GenomeConfig gcon) {
		econf = econ;
		gconf = gcon;
	}
	
	public void execute(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager =  new ExperimentManager(econf);
		float[][] densities = new float[regions.size()][manager.getSamples().size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		int SampleCount = 0;
		for(Sample sample : manager.getSamples()){
			double rpm = sample.getHitCount()/1000000; // Reads per million
			header.append(sample.getName());
			header.append("\t");
			int RegionCount = 0;
			for(Region r : regions){
				densities[RegionCount][SampleCount] = (float)(sample.countHits(r)/r.getWidth()*rpm);
				RegionCount++;
			}
			SampleCount++;
		}
		
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		
		StringBuilder densitiesSb = new StringBuilder();
		for(int r=0; r<regions.size(); r++){
			densitiesSb.append(regions.get(r).getLocationString());
			densitiesSb.append("\t");
			for(int d=0; d<manager.getSamples().size(); d++){
				densitiesSb.append(densities[r][d]);
				densitiesSb.append("\t");
			}
			densitiesSb.deleteCharAt(densitiesSb.length()-1);
			densitiesSb.append("\n");
		}
		
		System.out.println(densitiesSb.toString());
		manager.close();
		
		
	}
	
	//Mutators
	
	public void setRegions(List<Region> regs){regions=regs;}
	
	
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(),args);
		
		ArgParser ap = new ArgParser(args);
		ChipSeqExptRegionDensityProfiler profiler = new ChipSeqExptRegionDensityProfiler(econ,gcon);
		
		int win  = Args.parseInteger(args, "win", 1000);
		
		if(!ap.hasKey("peaks") && !ap.hasKey("regions")){
			System.err.println("Provide peaks or regions file!!!");
			System.exit(1);
		}else{
			if(ap.hasKey("peaks")){
				List<Point> points = RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), ap.getKeyValue("peaks"), -1);
				List<Region> regs = new ArrayList<Region>();
				for(Point p : points){
					regs.add(p.expand(win));
				}
				profiler.setRegions(regs);
			}else{
				List<Region> regs = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), ap.getKeyValue("regions"), -1);
				profiler.setRegions(regs);
			}
			
			profiler.execute();
			
		}
		
		
		
	}
	
	
	
	
	
	
	
	
	

}
