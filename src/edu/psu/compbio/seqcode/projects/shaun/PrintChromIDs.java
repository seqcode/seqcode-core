package edu.psu.compbio.seqcode.projects.shaun;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class PrintChromIDs {

	public static void main(String[] args) {
		
		try {
			Pair<Species, Genome> pair = Args.parseGenome(args);
	    	Species currorg = pair.car();
	    	Genome currgen = pair.cdr();
	    	
	    	for(String chr : currgen.getChromList()){
	    		System.out.println(currgen.getVersion()+"\t"+currgen.getDBID()+"\t"+chr+"\t"+currgen.getChromID(chr));
	    	}
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
