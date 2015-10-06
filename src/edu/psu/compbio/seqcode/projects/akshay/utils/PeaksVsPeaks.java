package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class PeaksVsPeaks {
	
	protected Genome gen;
	protected List<Point> peaksA;
	protected List<Region> regionsA;
	protected List<Point> peaksB;
	protected List<Region> regionsB;
	
	protected int overlapD; // min distance between the closest peaks to call them overlapping 
	protected int win; // for expanding a peak into a peak region	
	
	public PeaksVsPeaks(Genome g) {
		gen = g;
	}
	
	
	
	// Settors
	public void setPeaksA(List<Point> a){peaksA = a;}
	public void setPeaksB(List<Point> b){peaksB = b;}
	public void setRegsA(List<Region> a){regionsA = a;}
	public void setRegsB(List<Region>  b){regionsB = b;}
	public void setOverlapD(int o){overlapD =o;}
	public void setWin(int w){win = w;}
	
	
	public void printClosestPeaks(){
		// Need not hash by chrom. Usually lists are very small
		for(Point pa : peaksA){
			int mind=Integer.MAX_VALUE;
			Point nearestPeak=null;
			for(Point pb : peaksB){
				if(pb.getChrom().equals(pa.getChrom())){
					if(pa.distance(pb) < mind){
						mind = pa.distance(pb);
						nearestPeak = pb;
					}
				}
			}
			if(mind <= overlapD){
				System.out.println(pa.getLocationString()+"\t"+nearestPeak.getLocationString()+"\t"+Integer.toString(mind));
			}else{
				System.out.println(pa.getLocationString()+"\t-");
			}
			
		}
		
	}
	
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		
		String peaksAfile = Args.parseString(args, "peaksA", null);
		if(peaksAfile == null){
			System.err.println("Provide ChIP-Seq peaks file!!");
			return;
		}
		
		String peaksBfile = Args.parseString(args, "peaksB", null);
		if(peaksBfile == null){
			System.err.println("Provide ChIP-Seq peaks file!!");
			return;
		}
		
		int win = Args.parseInteger(args, "win", 150);
		
		List<Point> peaksA = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		List<Point> peaksB = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksBfile, win);
		
		List<Region> regionsA = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		List<Region> regionsB = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksBfile, win);
		
		
		PeaksVsPeaks analyzer = new PeaksVsPeaks(gconfig.getGenome());
		
		analyzer.setPeaksA(peaksA);
		analyzer.setPeaksB(peaksB);
		analyzer.setRegsA(regionsA);
		analyzer.setRegsB(regionsB);
		
		int overlapD = Args.parseInteger(args, "overlapD", 100);
		analyzer.setOverlapD(overlapD);
		
		
		analyzer.printClosestPeaks();
	}
	
	
	
	

}
