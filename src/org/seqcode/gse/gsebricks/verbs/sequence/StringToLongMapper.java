package org.seqcode.gse.gsebricks.verbs.sequence;

import org.seqcode.gse.gsebricks.verbs.Mapper;
import org.seqcode.gse.utils.sequence.SequenceUtils;

public class StringToLongMapper implements Mapper<String,Long> {

    public Long execute(String a) {
        return SequenceUtils.StringToLong(a);
    } 	   
}