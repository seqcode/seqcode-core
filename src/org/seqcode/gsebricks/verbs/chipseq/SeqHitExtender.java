package org.seqcode.gsebricks.verbs.chipseq;

import org.seqcode.data.seqdata.SeqHit;
import org.seqcode.gsebricks.verbs.Mapper;

public class SeqHitExtender implements Mapper<SeqHit,SeqHit> {
	
	private int extension;
	
	public SeqHitExtender(int ext) { 
		extension = ext;
	}

	public SeqHit execute(SeqHit a) {
        return a.extendHit(extension);
	}
}
