package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author mahony
 * Represents an entry in the core.readtype table.  
 * This is used as metadata for sequencing experiments.  
 * The name reflects the type of sequencing read, which is "SINGLE" or "PAIRED" for now
 *
 */
public class ReadType implements Comparable<ReadType>{
	private int dbid;
	private String name;
    
    /**
     * Creates a new <code>ReadType</code> object from
     * a ResultSet that is set to a row with two
     * values:
     *  an integer that is the database id
     *  a string that is the name
     */
    public ReadType(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
    }
	
	public String getName() { return name; }
	public int getDBID() { return dbid; }
    
    public String toString() { return name + " (#" + dbid + ")"; }
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ReadType)) { return false; }
		ReadType c = (ReadType)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(ReadType r) { return name.compareTo(r.name); }
}
