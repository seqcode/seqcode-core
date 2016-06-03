package org.seqcode.gse.gsebricks.verbs.chipseq;

import org.seqcode.gse.datasets.seqdata.SeqHit;
import org.seqcode.gse.gsebricks.verbs.Mapper;

public class SeqHitExtender implements Mapper<SeqHit,SeqHit> {
	
	private int extension;
	
	public SeqHitExtender(int ext) { 
		extension = ext;
	}

	public SeqHit execute(SeqHit a) {
        return a.extendHit(extension);
	}
}
