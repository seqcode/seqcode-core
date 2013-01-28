package edu.psu.compbio.seqcode.gse.ewok.verbs;

public class DefaultToStringMapper implements Mapper<Object, String> {

    public String execute(Object o) {
        return o.toString();
    }

}