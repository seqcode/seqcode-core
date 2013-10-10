package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

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
	
	public List<Region> parsePeaksToRegionsWithOutMerging(Collection<String> files, int win){
		List<Region> regs = new ArrayList<Region>();
		List<Point> points = new ArrayList<Point>();
		//Load the points
		for(String f : files){
			points.addAll(Utils.loadPointsFromFile(f, config.getGenome()));
		}
		for(Point p : points){
			Region r = p.expand(win/2);
			regs.add(r);
		}
		//Simple stats
		int totalLen = 0;
		for(Region r : regs)
			totalLen += r.getWidth();
		System.out.println(regs.size()+" regions with "+totalLen+"bp total length.");
				
		return regs;
	}
	
	/**
	 * Parse a set of region files
	 * @param files
	 * @param win
	 * @return
	 */
	public List<Region> parseRegions(Collection<String> files, int win){
		List<Region> regs = new ArrayList<Region>();
		
		for(String f : files){
			regs.addAll(Utils.loadRegionsFromFile(f, config.getGenome(), win));
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
	 * Count reads in the peakRegions for all given replicates
	 * @param peakRegions
	 */
	public void countSignal(List<Region> peakRegions){
		ExperimentSet eset = manager.getExperimentSet();
		if(eset.getConditions().size()==0){
			System.out.println("No experiments specified."); System.exit(1);
		}
		for(ExperimentCondition c : eset.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				double total = r.getSignal().getHitCount();
				System.out.println("Condition "+c.getName()+":\tRep "+r.getName());
				System.out.println("\tTotal Hit Count:\t"+total);
				
				for(Region pr : peakRegions){
					double currSig = r.getSignal().countHits(pr);
					System.out.println(pr.getLocationString()+"\t"+currSig);
				}
			}
		}
	}
	
	public void countSignalAcrossReplicates(List<Region> peakRegions){
		ExperimentSet eset = manager.getExperimentSet();
		if(eset.getConditions().size()==0){
			System.out.println("No experiments specified."); System.exit(1);
		}
		
		// print header line
		String header ="Point"+"\t";
		for(ExperimentCondition c :eset.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				header = header + r.getName()+"\t";
			}
		}
		System.out.println(header);
		
		for(Region pr : peakRegions){
			String out="";
			out = pr.getMidpoint().getChrom()+":"+pr.getMidpoint().getLocation();
			for(ExperimentCondition c : eset.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					double currSig = r.getSignal().countHits(pr);
					out=out+"\t"+Double.toString(currSig);
				}
			}
			System.out.println(out);
		}
	}
	
	/**
	 * Print quadrant values around stranded points
	 * @param motifhits
	 * @param quadrants
	 */
	public void printQuadValues(List<StrandedPoint> motifhits, Integer[] quadrants){
		ExperimentSet eset = manager.getExperimentSet();
		if(eset.getConditions().size()==0){
			System.out.println("No experiments specified."); System.exit(1);
		}
		
		int max=0; 
		for(Integer q : quadrants){
			if(Math.abs(q)>max)
				max = Math.abs(q);
		}
		
		//All loaded experiments are flattened
		for(StrandedPoint sp : motifhits){
			Region currReg = sp.expand(max);
			double[] pos_hits = new double[(max*2)+1];
			double[] neg_hits = new double[(max*2)+1];
			for(int h=0; h<pos_hits.length; h++){
				pos_hits[h]=0; neg_hits[h]=0;
			}
			
			for(ExperimentCondition c : eset.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					List<StrandedBaseCount> bc = r.getSignal().getUnstrandedBases(currReg);
					for(StrandedBaseCount b : bc){
						if(b.getStrand()=='+' && sp.getStrand()=='+')
							pos_hits[b.getCoordinate()-currReg.getStart()]+=b.getCount();
						else if(b.getStrand()=='-' && sp.getStrand()=='-')
							pos_hits[currReg.getEnd()-b.getCoordinate()]+=b.getCount();
						else if(b.getStrand()=='-' && sp.getStrand()=='+')
							neg_hits[currReg.getEnd()-b.getCoordinate()]+=b.getCount();
						else if(b.getStrand()=='+' && sp.getStrand()=='-')
							neg_hits[b.getCoordinate()-currReg.getStart()]+=b.getCount();
					}
				}
			}
			
			//Sum over quadrants
			int posQuad1=0, posQuad2=0, negQuad1=0, negQuad2=0;
			for(int i=quadrants[0]+max; i<quadrants[1]+max; i++){
				posQuad1+=pos_hits[i];
				negQuad1+=neg_hits[i];
			}
			for(int i=quadrants[1]+max; i<quadrants[2]+max; i++){
				posQuad2+=pos_hits[i];
				negQuad2+=neg_hits[i];
			}
			
			System.out.println(sp.getLocationString()+"\t"+posQuad1+"\t"+posQuad2+"\t"+negQuad2+"\t"+negQuad1);
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
					"\tOR\n" +
					"\t--regions <region file(s)>\n" +
					"\t--win <window size to take around each peak point>\n" +
					"\t--signalnoise [flag to print signal/noise stats]\n" +
					"\t--count [flag to count reads in signal regions]\n" +
					"\t--motifhits <motif hit file(s)>\n" +
					"\t--quadrants <intervals for upstream of site, comma-separated>\n" +
					"");
		}else{
			ExperimentManager manager = new ExperimentManager(config);
				
			SiteQuantifier quant = new SiteQuantifier(config, manager);
			
			//Load peak files if present and do signal/noise stats
			int win = Args.parseInteger(args,"win",200);
			if(ap.hasKey("peaks") || ap.hasKey("regions")){
				List<Region> peakRegions =null;
				if(ap.hasKey("peaks")){
					Collection<String> peakFiles = Args.parseStrings(args, "peaks");
					//peakRegions = quant.parsePeaksToRegions(peakFiles, win);
					peakRegions = quant.parsePeaksToRegionsWithOutMerging(peakFiles, win);
				}
				if(ap.hasKey("regions")){
					Collection<String> peakFiles = Args.parseStrings(args, "regions");
					peakRegions = quant.parseRegions(peakFiles, -1);
				}
				
				if(Args.parseFlags(args).contains("signalnoise")){
					quant.calcSigNoiseRatios(peakRegions);
				}
				if(Args.parseFlags(args).contains("count")){
					//quant.countSignal(peakRegions);
					quant.countSignalAcrossReplicates(peakRegions);
				}
			}
			
			
			//Load motif hit file
			if(ap.hasKey("motifhits")){
				String motiffile = Args.parseString(args, "motifhits", null);
				List<StrandedPoint> spoints = Utilities.loadStrandedPointsFromMotifFile(config.getGenome(), motiffile, win);
			
				if(ap.hasKey("quadrants")){
					String quads = Args.parseString(args, "quadrants", null);
					String[] quadsB = quads.split(",");
					Integer[] quadrants = new Integer[quadsB.length];
					for(int i=0; i<quadsB.length; i++)
						quadrants[i] = new Integer(quadsB[i]);
					
					quant.printQuadValues(spoints, quadrants);
				}
			}
			
			
			manager.close();
		}
	}

}
