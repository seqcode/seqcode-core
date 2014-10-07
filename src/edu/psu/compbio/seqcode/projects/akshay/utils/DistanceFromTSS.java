package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class DistanceFromTSS {
	
	Genome gen;
	List<StrandedPoint> refTSS = new ArrayList<StrandedPoint>();
	List<Point> locations = new ArrayList<Point>();
	
	public DistanceFromTSS(Genome g, String refPts, String locs) {
		this.gen = g;
		this.refTSS = RegionFileUtilities.loadStrandedPointFromRefTssFile(gen, refPts);
		this.locations = RegionFileUtilities.loadPeaksFromPeakFile(gen, locs, 0);
	}
	
	public void printDistances(){
		for(Point p : this.locations){
			int dis = Integer.MAX_VALUE;
			for(StrandedPoint sp : this.refTSS){
				if(sp.getChrom().equals(p.getChrom())){
					
					int tempd = sp.distance(p);
					if(tempd < dis){
						dis = tempd;
					}
				}
			}
			System.out.println(p.getLocationString()+"\t"+Integer.toString(dis));
		}
	}
	
	
	
	public static void main(String[] args) throws NotFoundException{
		ArgParser ap = new ArgParser(args);
		Genome g = Args.parseGenome(args).cdr();
		String tss = ap.getKeyValue("ref");
		String locs = ap.getKeyValue("peaks");
		
		DistanceFromTSS analyzer =  new DistanceFromTSS(g,tss,locs);
		System.out.println("Distance");
		analyzer.printDistances();
		
	}

}
