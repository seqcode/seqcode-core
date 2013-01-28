package edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;

public class ChipSeqHitExtender implements Mapper<ChipSeqHit,ChipSeqHit> {
	
	private int extension;
	
	public ChipSeqHitExtender(int ext) { 
		extension = ext;
	}

	public ChipSeqHit execute(ChipSeqHit a) {
        return a.extendHit(extension);
	}
}
