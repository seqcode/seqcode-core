package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

public class DistanceFromTSS {
	
	Genome gen;
	List<StrandedPoint> refTSS = new ArrayList<StrandedPoint>();
	List<Point> locations = new ArrayList<Point>();
	
	public DistanceFromTSS(Genome g, String refPts, String locs) {
		this.gen = g;
		this.refTSS = Utilities.loadStrandedPointFromRefTssFile(gen, refPts);
		this.locations = Utilities.loadPeaksFromPeakFile(gen, locs, 0);
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
			System.out.println(p.getLocation()+"\t"+Integer.toString(dis));
		}
	}
	
	
	
	public static void main(String[] args) throws NotFoundException{
		ArgParser ap = new ArgParser(args);
		Genome g = Args.parseGenome(args).cdr();
		String tss = ap.getKeyValue("ref");
		String locs = ap.getKeyValue("peaks");
		
		DistanceFromTSS analyzer =  new DistanceFromTSS(g,tss,locs);
		analyzer.printDistances();
		
	}

}
