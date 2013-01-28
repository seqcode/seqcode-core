package edu.psu.compbio.seqcode.projects.shaun;

import java.util.*;
import java.util.regex.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ExpanderIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Filter;
import edu.psu.compbio.seqcode.gse.ewok.verbs.TiledRegionGenerator;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class WellTiledRegionParser {
	
	public static void main(String[] args) { 
		Genome g = null;
		try {
			g = Organism.findGenome(args[0]);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		WellTiledRegionParser wtp = new WellTiledRegionParser(g, args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]));
		try {
			wtp.save(new File(args[4]));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private int spacing;
	private Genome genome;
	private Set<Region> regions;

	public WellTiledRegionParser() throws IOException { 
		File f = new File("well_tiled_regions.txt");
		parse(f);
		spacing = 10000;
	}
	
	public WellTiledRegionParser(File f, int sp) throws IOException {
		parse(f);
		spacing = sp;
	}
	
	public WellTiledRegionParser(Genome g, String arrayName, int sp, int mincount) { 
		genome = g;
		regions = new HashSet<Region>();
		this.spacing = 10000; 
		try {
			TiledRegionGenerator<NamedRegion> tiler = new TiledRegionGenerator<NamedRegion>(arrayName, sp, mincount);
			ChromRegionIterator chroms = new ChromRegionIterator(genome);
			Iterator<Region> tileditr = new ExpanderIterator<NamedRegion,Region>(tiler, chroms);
			
			while(tileditr.hasNext()) { 
				regions.add(tileditr.next());
			}
			
		} catch (NotFoundException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(arrayName);
		}
	}
	
	public void save(File f) throws IOException { 
		PrintStream ps = new PrintStream(new FileOutputStream(f));
		for(Region r : regions) { 
			ps.println(String.format("%s\t%d-%d", r.getChrom(), r.getStart(), r.getEnd()));
		}
		ps.close();
	}
		
	public void parse(File f) throws IOException { 
		BufferedReader br = new BufferedReader(new FileReader(f));
		Pattern p = Pattern.compile("([^\\s]+)\\s+(\\d+)[-\\s]+(\\d+)");
		try {
			genome = Organism.findGenome("mm8");
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		regions = new HashSet<Region>();
		
		String line;
		while((line = br.readLine()) != null) { 
			if((line=line.trim()).length() > 0) { 
				Matcher m = p.matcher(line);
				if(m.matches()) { 
					String chrom = m.group(1);
					int start = Integer.parseInt(m.group(2));
					int end = Integer.parseInt(m.group(3));
					Region r = new Region(genome, chrom, start, end);
					regions.add(r);
				}
			}
		}
		
		br.close();
	}
	
	public boolean isWellTiled(Region r) { 
		for(Region t : regions) { 
			if(t.contains(r)) { 
				return true;
			}
		}
		return false;
	}
	
	public boolean isWellTiled(Gene g) { 
		int tss = g.getTSS();
		Region r = new Region(genome, g.getChrom(), tss-spacing, tss+spacing);
		return isWellTiled(r);
	}

	public Set<Region> getRegions() {
		return regions;
	}
    
    public Filter<Region,Region> createFilter() { 
        return new Filter<Region,Region>() {
            public Region execute(Region a) {
                return isWellTiled(a) ? a : null;
            } 
        };
    }
}
