package edu.psu.compbio.seqcode.gse.gsebricks.verbs.chipseq;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Mapper;

public class SeqHitExtender implements Mapper<SeqHit,SeqHit> {
	
	private int extension;
	
	public SeqHitExtender(int ext) { 
		extension = ext;
	}

	public SeqHit execute(SeqHit a) {
        return a.extendHit(extension);
	}
}
