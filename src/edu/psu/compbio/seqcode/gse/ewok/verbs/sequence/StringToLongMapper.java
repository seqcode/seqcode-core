package edu.psu.compbio.seqcode.gse.ewok.verbs.sequence;

import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class StringToLongMapper implements Mapper<String,Long> {

    public Long execute(String a) {
        return SequenceUtils.StringToLong(a);
    } 	   
}