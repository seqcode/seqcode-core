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
	List<Point> refTSS = new ArrayList<Point>();
	List<Point> locations = new ArrayList<Point>();
	
	public DistanceFromTSS(Genome g, String refPts, String locs) {
		this.gen = g;
		this.refTSS = Utilities.loadPeaksFromPeakFile(gen, refPts, 0);
		this.locations = Utilities.loadPeaksFromPeakFile(gen, locs, 0);
	}
	
	public void printDistances(){
		for(Point p : this.locations){
			int dis = Integer.MAX_VALUE;
			for(Point sp : this.refTSS){
				int tempd=0;
				if(sp.getChrom().equals(p.getChrom()))
					tempd = sp.distance(p);
					if(tempd < dis){
						dis = tempd;
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
