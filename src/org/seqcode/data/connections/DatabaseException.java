package org.seqcode.data.connections;

public class DatabaseException extends RuntimeException {

    public DatabaseException(String s) {
        super(s);
    }
    public DatabaseException(String s, Exception e) {
        super(s,e);
    }

}
