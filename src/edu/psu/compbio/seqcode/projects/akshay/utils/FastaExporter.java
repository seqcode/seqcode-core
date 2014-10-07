package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.EventMetaMaker;

public class FastaExporter {
	
	
	public List<Region> regs;
	public Genome gen;
	public int win;
	public SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();;
	public List<String> fasta = new ArrayList<String>();
	
	
	public FastaExporter(int win, List<Point> points, Genome gen, boolean cache, String seqPathName) {
		regs = new ArrayList<Region>();
		this.win = win;
		this.gen = gen;
		for(Point p: points){
			Region reg = p.expand(win/2);
			regs.add(reg);
		}
		if(cache){
			seqgen.useCache(cache);
			seqgen.setGenomePath(seqPathName);
		}
		
		for(Region r : regs){
			String fa = seqgen.execute(r);
			fasta.add(fa);
		}
	}
	
	public void printFasta(){
		for(int i=0; i< regs.size(); i++){
			System.out.println(">"+regs.get(i).getLocationString());
			System.out.println(fasta.get(i));
		}
	}
	
	
	public static void main(String[] args){
		try{
			ArgParser ap =  new ArgParser(args);
			int window = Args.parseInteger(args, "win", 200);
			boolean cache = ap.hasKey("cache");
			String seqPathName = "";
			if(cache){
				seqPathName = Args.parseString(args, "seq", "");
			}
			String points = ap.getKeyValue("locations");
			
			Genome gen;
			if(ap.hasKey("species") || ap.hasKey("genome") || ap.hasKey("gen")){
				Pair<Organism, Genome> orggen = Args.parseGenome(args);
				gen = orggen.cdr();
			} 
			else{
				if(ap.hasKey("geneinfo") || ap.hasKey("g")){
					String infofilename = ap.hasKey("geneinfo") ? ap.getKeyValue("geneinfo") : ap.getKeyValue("g"); 
					gen = new Genome("Genome",new File(infofilename),true);
				}
				else{
					gen = null;
				}
			}
			List<Point> search_locs = RegionFileUtilities.loadPeaksFromPeakFile(gen, points, window);
			FastaExporter exporter = new FastaExporter(window, search_locs, gen, cache, seqPathName);
			exporter.printFasta();
		}catch(NotFoundException e){
			e.printStackTrace();
		} 
	}

}
