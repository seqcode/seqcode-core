package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.ArrayList;
import java.util.Arrays;
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
	
	protected List<List<Point>> peaksBSet;
	protected List<List<Region>> regionsBSet;
	
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
	public void setPeaksBSet(List<List<Point>> pset){peaksBSet = pset;}
	public void setRegionsBSet(List<List<Region>> rset){regionsBSet = rset;}
	
	
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
	
	public void printPeakSets(){
		for(Point pa: peaksA){
			StringBuilder sb= new StringBuilder();
			sb.append(pa.getLocationString());sb.append("\t");
			Point[] nearestPeaks = new Point[peaksBSet.size()];
			int[] nearestDis = new int[peaksBSet.size()];
			Arrays.fill(nearestDis, Integer.MAX_VALUE);
			int setId=0;
			for(List<Point> bset: peaksBSet){
				for(Point pb : bset){
					if(pa.distance(pb) < nearestDis[setId]){
						nearestDis[setId] = pa.distance(pb);
						nearestPeaks[setId] =pb;
					}
				}
				setId++;
			}
			for(int i=0; i<nearestDis.length;i++){
				if(nearestDis[i]<=overlapD){
					sb.append(nearestPeaks[i].getLocationString());sb.append("\t");
				}else{
					sb.append(nearestPeaks[i].getLocationString());sb.append("\t");
				}
			}
			sb.deleteCharAt(sb.length()-1);
			System.out.println(sb.toString());
		}
	}
	
	
	
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		
		int win = Args.parseInteger(args, "win", 150);
		
		PeaksVsPeaks analyzer = new PeaksVsPeaks(gconfig.getGenome());
		
		String peaksAfile = Args.parseString(args, "peaksA", null);
		if(peaksAfile == null){
			System.err.println("Provide ChIP-Seq peaks file A!!");
			return;
		}
		List<Point> peaksA = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		List<Region> regionsA = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		analyzer.setPeaksA(peaksA);
		analyzer.setRegsA(regionsA);
		
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksSetB")){
			System.err.println("Provide ChIP-Seq peaks file to compare with the A list!!");
			return;
		}
		
		String peaksBfile = null;
		if(ap.hasKey("peaksB")){
			peaksBfile = Args.parseString(args, "peaksB", null);
			List<Point> peaksB = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksBfile, win);
			List<Region> regionsB = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksBfile, win);
			analyzer.setPeaksB(peaksB);
			analyzer.setRegsB(regionsB);
		}
		
		String[] peaksBfileList = null;
		if(ap.hasKey("peaksBset")){
			peaksBfileList =  ap.getKeyValue("peaksBset").split(";");
			List<List<Point>> bPeaksSet = new ArrayList<List<Point>>();
			List<List<Region>> bRegionSet = new ArrayList<List<Region>>();
			for(int s=0; s<peaksBfileList.length; s++){
				List<Point> currPoints = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksBfileList[s], win);
				bPeaksSet.add(currPoints);
				List<Region> currRegions = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksBfileList[s], win);
				bRegionSet.add(currRegions);
			}
			analyzer.setPeaksBSet(bPeaksSet);
			analyzer.setRegionsBSet(bRegionSet);
		}
		
		
		int overlapD = Args.parseInteger(args, "overlapD", 100);
		analyzer.setOverlapD(overlapD);
		
		if(ap.hasKey("peaksB"))
			analyzer.printClosestPeaks();
		if(ap.hasKey("peaksBset"))
			analyzer.printPeakSets();
	}
	
	
	
	

}
