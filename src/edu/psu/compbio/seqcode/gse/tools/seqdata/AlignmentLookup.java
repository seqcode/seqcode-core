package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class AlignmentLookup {

	public static void main(String[] args) {
		if(Args.parseArgs(args).contains("expt") && Args.parseArgs(args).contains("species")){
			SeqDataLoader loader=null;
			try {		
				Pair<Organism,Genome> pair = Args.parseGenome(args);
				Genome gen = pair.cdr();
				List<SeqLocator> expts = Args.parseSeqExpt(args,"expt");
				if(expts.size()>1){
					System.err.println("Specify just one alignment.");
					System.exit(1);
				}
				SeqLocator expt = expts.get(0);
				loader = new SeqDataLoader(false);
			
				if(loader.loadAlignments(expt, gen).size()>0)
					System.out.println("FOUND");
				else
					System.out.println("NOTFOUND");
		
				gen.close();
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
