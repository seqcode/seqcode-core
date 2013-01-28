package edu.psu.compbio.seqcode.gse.utils.database;

public class DatabaseException extends RuntimeException {

    public DatabaseException(String s) {
        super(s);
    }
    public DatabaseException(String s, Exception e) {
        super(s,e);
    }

}
