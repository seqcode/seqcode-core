package edu.psu.compbio.seqcode.gse.utils.database;

public class UnknownRoleException extends RuntimeException {
    public UnknownRoleException(String s) {
        super(s);
    }
    public UnknownRoleException(String s, Exception e) {
        super(s,e);
    }

}
