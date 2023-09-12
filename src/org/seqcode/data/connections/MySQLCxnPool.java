package org.seqcode.data.connections;

import java.util.*;

//         connectString = "jdbc:mysql://opteron.csail.mit.edu/" + dbname + "?user=" + username +
//-             "&password=" + password;



public class MySQLCxnPool extends CxnPool {

    public MySQLCxnPool(Properties props) {
        super(props);
        try {
            Class.forName("com.mysql.cj.jdbc.Driver").newInstance();
        } catch (Exception ex) {
            System.err.println("Couldn't load MySQL JDBC driver");
        }
    }
    public int getType() {return DatabaseFactory.MYSQL;}
}
        
