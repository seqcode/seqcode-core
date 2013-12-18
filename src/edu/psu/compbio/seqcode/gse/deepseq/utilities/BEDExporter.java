package edu.psu.compbio.seqcode.gse.deepseq.utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.apache.log4j.Logger;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.discovery.SingleConditionFeatureFinder;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * Outputs a BED format file for a deep-seq experiment.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class BEDExporter {
	private Organism org;
	private Genome gen;
	protected DeepSeqExpt expt;
	protected int [] stackedHitCountsPos;
	protected int [] stackedHitCountsNeg;
	private String outName="out";
	private boolean dbconnected=false;
	private static final Logger logger = Logger.getLogger(SingleConditionFeatureFinder.class);
	
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		
		BEDExporter bed = new BEDExporter(args);
		bed.execute();
		
		System.exit(0);
	}
	
	
	public BEDExporter(String [] args) {
		if(args.length==0){
			System.err.println("BEDExporter usage:\n" +
					"\t--species <organism;genome>\n" +
					"\t--(rdb)expt <experiment names>\n" +
					"\t--out <output file name>");
			System.exit(1);
		}
		ArgParser ap = new ArgParser(args);
		try {
			if(ap.hasKey("species")){
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				if(pair != null){
					gen = pair.cdr();
					dbconnected=true;
				}
			}else{
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("geninfo") || ap.hasKey("g")){
					String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
					gen = new Genome("Genome", new File(fName), true);
				}else{
				    gen = null;
				}
			}
		}catch (NotFoundException e) {
			e.printStackTrace();
		}
		
		outName = Args.parseString(args,"out",outName);
	    		
	    // Load the experiments
	    List<SeqLocator> dbexpts = Args.parseSeqExpt(args, "dbexpt");
	    List<SeqLocator> rdbexpts = Args.parseSeqExpt(args,"rdbexpt");
	    List<File> expts = Args.parseFileHandles(args, "expt");
	    boolean nonUnique = ap.hasKey("nonunique") ? true : false;
	    String fileFormat = Args.parseString(args, "format", "ELAND");
	    if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
	       	expt= new DeepSeqExpt(gen, expts, nonUnique, fileFormat, 1);
	    }else if (dbexpts.size() > 0 && expts.size() == 0) {
	    	expt = new DeepSeqExpt(gen, dbexpts, "db", 1);
	    	dbconnected = true;
	    }else if (rdbexpts.size()>0 && expts.size() == 0){
	    	expt = new DeepSeqExpt(gen, rdbexpts, "readdb", -1);
	    	dbconnected=true;
	    }else {
	      logger.error("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");
	      System.exit(1);
	    }
	    logger.info("Expt hit count: " + (int) expt.getHitCount() + ", weight: " + (int) expt.getWeightTotal());
	}
	
	public void execute(){
		try {
			FileWriter fw = new FileWriter(outName);
						       			
			ChromRegionIterator chroms = new ChromRegionIterator(gen);
			while(chroms.hasNext()){
				NamedRegion currentRegion = chroms.next();
				
				//Split the job up into chunks of 100Mbp
				for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=100000000){
					int y = x+100000000; 
					if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
					Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
					
					String BEDhits = expt.getBED_StrandedReads(currSubRegion, '+', 2.0);
                    fw.write(BEDhits);
                    BEDhits = expt.getBED_StrandedReads(currSubRegion, '-', 2.0);
                    fw.write(BEDhits);
				}
			}
			
			fw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	//keep the number in bounds
	protected final double inBounds(double x, double min, double max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
}
