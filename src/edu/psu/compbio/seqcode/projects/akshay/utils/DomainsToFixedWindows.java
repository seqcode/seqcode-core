package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class DomainsToFixedWindows {
	
	private List<Region> domains;
	private int win;
	private boolean subset;
	private GenomeConfig gconf;
	
	
	public DomainsToFixedWindows(GenomeConfig gc) {
		gconf = gc;
	}
	
	//Mutators
	public void setDomains(List<Region> doms){domains = doms;}
	public void setWin(int w){win=w;}
	public void selectSubset(boolean sub){subset = sub;}
	
	
	public void execute(){
		for(Region domain: domains){
			if(domain.getWidth() < win*2){
				System.out.println(domain.getMidpoint());
			}else{
				int domain_start = domain.getStart();
				int currLocation = domain_start+win/2;
				while(currLocation+win/2<=domain.getEnd()){
					Point p = new Point(gconf.getGenome(),domain.getChrom(),currLocation);
					System.out.println(p.getLocationString());
					currLocation = currLocation+win;
				}
			}
		}
	}
	
	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		GenomeConfig gc = new GenomeConfig(args);
		int win = Args.parseInteger(args, "win", 150);
		boolean sub = ap.hasKey("subset");
		String domains_filename = ap.getKeyValue("domains");
		List<Region> doms = RegionFileUtilities.loadRegionsFromPeakFile(gc.getGenome(), domains_filename, -1);
		
		DomainsToFixedWindows runner = new DomainsToFixedWindows(gc);
		runner.setDomains(doms);
		runner.setWin(win);
		runner.selectSubset(sub);
		
		runner.execute();
		
		
		
	}
	
	

}
