package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class PairEndStatSandbox {
	
	// Current this class uses the old seqexpander methods instead of the new hitloaders. ..
	// Change to the new hitloader methods as soon as possible
	
	private SeqLocator locator;
	
	private LinkedList<SeqAlignment> alignments;
	private Genome lastGenome;
	private SeqDataLoader loader;
	private boolean closeLoader;
	private Genome gen;
	
	public PairEndStatSandbox(SeqLocator expt, Genome g) {
		this.locator = expt;
		
		try {
			this.getAligns(g);
		} catch (SQLException e) {
			e.printStackTrace();
		}	
	}
	
	//private Iterator<SeqHitPair> loadHitsbyChrom(){
	//	gen.getChromList()
	//	loader.loadPairsByChrom(a, chromid)
	//}
	
	
	
	
	private void getAligns(Genome genome) throws SQLException {
		if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<SeqAlignment>();
        lastGenome = genome;
        try {
            alignments.addAll(locator.loadAlignments(loader, genome));
        } catch (SQLException e) {
            e.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
	
	public void close() {
        if (closeLoader) {
            loader.close();
        }
        loader = null;
        if (alignments != null) {
            alignments.clear();
        }
    }


    public boolean isClosed() {
        return loader == null;
    }
	
	
	public static void main(String[] args){
		SeqLocator expt = Args.parseSeqExpt(args, "expt").get(0);
		
	}
}
