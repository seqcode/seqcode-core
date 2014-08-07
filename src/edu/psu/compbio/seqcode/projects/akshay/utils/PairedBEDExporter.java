package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.projects.readdb.ClientException;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class PairedBEDExporter {

	
	private Genome gen;
	
	private List<SeqLocator> locators;
	
	private String outName="out";
	private DeepSeqExpt fetcher;
	
	public PairedBEDExporter(List<SeqLocator> expt, Genome g, String o) throws SQLException, IOException, NotFoundException {
		this.gen = g;
		this.locators = expt;
		this.fetcher = new DeepSeqExpt(gen,this.locators,"readdb",-1,true);
		this.outName =o;
	}

	
	public void execute() throws IOException{
		
		FileWriter fw = new FileWriter(outName);
		ChromRegionIterator chroms = new ChromRegionIterator(gen);
			
		while(chroms.hasNext()){
			NamedRegion currentRegion =  chroms.next();
				
			for(int x= currentRegion.getStart(); x<= currentRegion.getEnd(); x+=100000000){
				int y= x+100000000;
				if(y>currentRegion.getEnd()){y = currentRegion.getEnd();}
				Region currSubRegion = new Region(gen, currentRegion.getChrom(),x,y);
					
				String BEDhits = this.getMinPointBed(currSubRegion);
				fw.write(BEDhits);
			}
			
		}
		
		fw.close();
	}
	
	private List<PairedHit> getPairs(Region r){
		List<PairedHit> ret = this.fetcher.loadPairsAsPairs(r);
		return ret;
	}
	
	public String getMinPointBed(Region r){
		String head = "chr"+r.getChrom()+"\t";
		String tail="\tU\t0\t"+"+";
		StringBuilder sb = new StringBuilder();
		for(PairedHit ph : this.getPairs(r)){
			int midp=0;
			if(ph.leftChrom == ph.rightChrom){
				midp = ph.lesserPos() + (ph.greaterPos()-ph.lesserPos())/2;
				
				sb.append(head).append(midp-1).append("\t").append(midp).append(tail).append("\n");
			}
		}
		return sb.toString();
	}
	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		Genome g=null;
		try {
			if(ap.hasKey("species")){
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				if(pair != null){
					g = pair.cdr();
					
				}
			}else{
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("geninfo") || ap.hasKey("g")){
					String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
					g = new Genome("Genome", new File(fName), true);
				}else{
				    g = null;
				}
			}
		}catch (NotFoundException e) {
			e.printStackTrace();
		}
		
		
		
		List<SeqLocator> rdbexpts = Args.parseSeqExpt(args, "rdbexpt");
		String out = ap.getKeyValue("out");
		try {
			PairedBEDExporter exporter = new PairedBEDExporter(rdbexpts,g,out);
			exporter.execute();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
}
