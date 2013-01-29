package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

/**
 * SiteQuantifier: 
 * 		- Count reads in signal regions vs background
 * 		- Split counts into windows around a motif hit
 *  
 * @author mahony
 *
 */
public class SiteQuantifier {
	
	protected ExperimentManager manager;
	protected Config config;
		
	public SiteQuantifier(Config con, ExperimentManager man){
		config = con;
		manager = man;
	}
	
	/**
	 * Parse a set of peak files, merging all points into expanded windows
	 * @param files
	 * @param win
	 * @return
	 */
	public List<Region> parsePeaksToRegions(Collection<String> files, int win){
		List<Region> regs = new ArrayList<Region>();
		List<Point> points = new ArrayList<Point>();
		//Load the points
		for(String f : files){
			points.addAll(Utils.loadPointsFromFile(f, config.getGenome()));
		}
		//Sort the points
		Collections.sort(points);
		//Convert points to regions
		Region lastr = null;
		for(Point p : points){
			//expand 
			Region r = p.expand(win/2);
			//Overlap check
			if(lastr != null && lastr.overlaps(r)){
				//Combine
				lastr = lastr.combine(r);
			}else{
				//Add to list
				regs.add(r);
				lastr = r;
			}
		}
		//Simple stats
		int totalLen = 0;
		for(Region r : regs)
			totalLen += r.getWidth();
		System.out.println(regs.size()+" regions with "+totalLen+"bp total length.");
		
		return regs;
	}
	
	/**
	 * Count reads in the peakRegions for all given replicates
	 * @param peakRegions
	 */
	public void calcSigNoiseRatios(List<Region> peakRegions){
		ExperimentSet eset = manager.getExperimentSet();
		if(eset.getConditions().size()==0){
			System.out.println("No experiments specified."); System.exit(1);
		}
		for(ExperimentCondition c : eset.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				double total = r.getSignal().getHitCount();
				double sig = 0;
				System.out.println("Condition "+c.getName()+":\tRep "+r.getName());
				System.out.println("\tTotal Hit Count:\t"+total);
				
				for(Region pr : peakRegions)
					sig+= r.getSignal().countHits(pr);
				
				System.out.println("\tSignal proportion = "+sig/total);
			}
		}
		
	}
	
	
	
	/**
	 * Main driver method for SiteQuantifier
	 * @param args
	 */
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("SiteQuantifier:");
			System.err.println("\tArgs:\n" +
					"\t--peaks <peaks file(s)>\n" +
					"\t--win <window size to take around each peak point>\n" +
					"\t--signalnoise [flag to print signal/noise stats]\n" +
					"");
		}else{
			ExperimentManager manager = new ExperimentManager(config);
				
			SiteQuantifier quant = new SiteQuantifier(config, manager);
			
			//Load peak files if present
			int win = Args.parseInteger(args,"win",200);
			if(ap.hasKey("peaks")){
				Collection<String> peakFiles = Args.parseStrings(args, "peaks");
				List<Region> peakRegions = quant.parsePeaksToRegions(peakFiles, win);
				
				if(Args.parseFlags(args).contains("signalnoise")){
					quant.calcSigNoiseRatios(peakRegions);
				}
			}
			
			
			
			
			manager.close();
		}
	}

}
