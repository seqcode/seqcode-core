package org.seqcode.gse.tools.seqdata;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.gse.datasets.seqdata.SeqDataLoader;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;


public class AlignmentLookup {

	public static void main(String[] args) {
		if(Args.parseArgs(args).contains("expt") && Args.parseArgs(args).contains("species")){
			SeqDataLoader loader=null;
			try {		
				Pair<Species,Genome> pair = Args.parseGenome(args);
				Genome gen = pair.cdr();
				List<SeqLocator> expts = Args.parseSeqExpt(args,"expt");
				if(expts.size()>1){
					System.err.println("Specify just one alignment.");
					System.exit(1);
				}
				SeqLocator expt = expts.get(0);
				loader = new SeqDataLoader(false, true);
			
				boolean erafound=false;
				boolean erfound=false;
				try{
					if(loader.loadAlignments(expt, gen).size()>0)
						erafound = true;
				} catch (NotFoundException e) {}
				
				try{
					if(loader.loadExperiment(expt.getExptName(), expt.getReplicateString())!=null)
						erfound=true;
				} catch (NotFoundException e) {}
				
				if(erafound)
					System.out.println("EXPTREPALIGNFOUND");
				else if(erfound)
					System.out.println("EXPTREPFOUND");
				else
					System.out.println("NOTFOUND");
			} catch (NotFoundException e) {
				System.out.println("NOTFOUND");
			} catch (SQLException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if(loader!=null)
					loader.close();
			}

    	}else{
    		System.err.println("AlignmentLookup:\n" +
    				"\t--species <species;genome>\n" +
    				"\t--expt <expt;rep;align>\n");
    	}
	}
}
