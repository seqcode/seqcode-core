package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;

public class PairEndStatSandbox {
	
	// Current this class uses the old seqexpander methods instead of the new hitloaders. ..
	// Change to the new hitloader methods as soon as possible
	
	SeqLocator locator;
	SeqExpander expander;
	
	
	public PairEndStatSandbox(SeqLocator expt) {
		this.locator = expt;
		try {
			this.expander = new SeqExpander(expt);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}
	
	
	
	
	
	public static void main(String[] args){
		SeqLocator expt = Args.parseSeqExpt(args, "expt").get(0);
		
	}
}
