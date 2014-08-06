package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.File;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

public class DnaMethOverlapStat {
	
	private List<Point> pts;
	private int winSize;
	private List<Region> regs;
	private Genome gen;
	private RRbsBEDLoader loader;
	
	
	
	public DnaMethOverlapStat(Genome g, String peaks_file, int win, String MethBED_File) {
		this.gen = g;
		this.winSize = win;
		this.pts = Utilities.loadPeaksFromPeakFile(g, peaks_file, this.winSize);
		this.regs = Utilities.loadRegionsFromPeakFile(gen, peaks_file, winSize);
		
		File Bed_file = new File(MethBED_File);
		this.loader = new RRbsBEDLoader(this.gen,Bed_file);
		loader.sourceReads();
	}
	
	
	public void printMethPercs(){
		for(Region r : this.regs){
			double perc = loader.getMethPerc(r).car();
			int num_sites = loader.getMethPerc(r).cdr();
			if(num_sites > 0){
				System.out.println(r.getLocationString()+"\t"+Double.toString(perc)+"\t"+num_sites);
			}
		}
	}
	
	
	public static void main(String[] args) throws NotFoundException{
		ArgParser ap = new ArgParser(args);
		
		String locations_file = ap.getKeyValue("points");
		String MethBed_file = ap.getKeyValue("methBed");
		
		Genome gen = null;
		if(ap.hasKey("species")){
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			if(pair != null){
				
				gen = pair.cdr();
			}
		}else{
			if(ap.hasKey("geninfo") || ap.hasKey("g")){
				//Make fake genome... chr lengths provided
				String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
				gen = new Genome("Genome", new File(fName), true);
			}else{
				gen = null;
			}
		}
		
		int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue() : 200;
		
		DnaMethOverlapStat analyzer = new DnaMethOverlapStat(gen, locations_file, win, MethBed_file);
		analyzer.printMethPercs();
		
	}
	
	

}
