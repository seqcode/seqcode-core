package edu.psu.compbio.seqcode.gse.datasets.core;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * SeqDataUser: represents an entry in the core.seqdatauser table.
 * Used as a convenient list of current users that have permission to view seqdata experiments.
 * Some users are administrators (admin=1).
 * The user "public" is treated as a special case in some database-interfacing methods. 
 * @author mahony
 *
 */
public class SeqDataUser implements Comparable<SeqDataUser>{
	private int dbid;
	private String name;
	private int admin;
    
    /**
     * Creates a new <code>SeqDataUser</code> object from
     * a ResultSet that is set to a row with three
     * values:
     *  an integer that is the database id
     *  a string that is the name
     *  an integer that is the admin flag
     */
    public SeqDataUser(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
        admin = rs.getInt(3);
    }
	
	public String getName() { return name; }
	public int getDBID() { return dbid; }
	public boolean isAdmin(){return admin==1;}
    
    public String toString() { return name; }
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof SeqDataUser)) { return false; }
		SeqDataUser c = (SeqDataUser)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(SeqDataUser r) { return name.compareTo(r.name); }
}
