package org.seqcode.projects.akshay.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;
import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.seqdata.SeqAlignment;
import org.seqcode.gse.datasets.seqdata.SeqDataLoader;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gse.projects.gps.DeepSeqExpt;
import org.seqcode.gse.projects.gps.discovery.SingleConditionFeatureFinder;
import org.seqcode.gse.projects.readdb.ClientException;
import org.seqcode.gse.projects.readdb.PairedHit;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;


public class PairedBEDExporter {

	
	private Genome gen;
	
	private List<SeqLocator> locators;
	
	private String outName="out";
	private DeepSeqExpt fetcher;
	
	private int minFragLen;
	private int maxFragLen;
	
	private boolean clientConnected = false;
	
	private static final Logger logger = Logger.getLogger(SingleConditionFeatureFinder.class);
	
	public PairedBEDExporter(List<SeqLocator> expt, Genome g, String o,int minFrag, int maxFrag) throws SQLException, IOException, NotFoundException {
		this.gen = g;
		this.locators = expt;
		this.fetcher = new DeepSeqExpt(gen,this.locators,"readdb",-1,false,true);
		this.clientConnected = true;
		this.logger.info("Expt hit count: " + (int) this.fetcher.getHitCount()+ ", weight: " + (int) this.fetcher.getWeightTotal());
		this.outName =o;
		this.maxFragLen =  maxFrag;
		this.minFragLen = minFrag;
	}
	
	public void close(){
		if(this.clientConnected){
			this.fetcher.closeLoaders();
		}
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
		this.close();
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
			int fraglen = ph.greaterPos() - ph.lesserPos(); 
			if(ph.leftChrom == ph.rightChrom && fraglen>=this.minFragLen && fraglen <= this.maxFragLen ){
				midp = ph.lesserPos() + (ph.greaterPos()-ph.lesserPos())/2;
				
				sb.append(head).append(midp-1).append("\t").append(midp).append(tail).append("\n");
			}
		}
		return sb.toString();
	}
	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		Genome g=null;
		int minFrag = 0;
		int maxFrag = 0;
		try {
			if(ap.hasKey("species")){
				Pair<Species, Genome> pair = Args.parseGenome(args);
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
			minFrag = Args.parseInteger(args, "minFragLen", 140);
			maxFrag = Args.parseInteger(args, "maxFragLen", 200);
			
		}catch (NotFoundException e) {
			e.printStackTrace();
		}
		
		List<SeqLocator>  rdbexpts = new ArrayList<SeqLocator>();
		if(ap.hasKey("reps")){
			String repstring = Args.parseString(args, "reps", "");
			Collection<String> reps = new ArrayList<String>();
			String[] pieces = repstring.split(";");
			for(int s=0; s<pieces.length; s++){
				reps.add(pieces[s]);
			}
			String exptname = Args.parseString(args, "rdbexpt", "");
			String[] exptPieces = exptname.split(";");
			rdbexpts.add(new SeqLocator(exptPieces[0],reps,exptPieces[1])); 
		}else{
			rdbexpts = Args.parseSeqExpt(args, "rdbexpt");
		}
		String out = ap.getKeyValue("out");
		try {
			PairedBEDExporter exporter = new PairedBEDExporter(rdbexpts,g,out,minFrag,maxFrag);
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
