package edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;

public class ChipSeqHitExtender implements Mapper<SeqHit,SeqHit> {
	
	private int extension;
	
	public ChipSeqHitExtender(int ext) { 
		extension = ext;
	}

	public SeqHit execute(SeqHit a) {
        return a.extendHit(extension);
	}
}
