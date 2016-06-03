package org.seqcode.projects.shaun;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;

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
