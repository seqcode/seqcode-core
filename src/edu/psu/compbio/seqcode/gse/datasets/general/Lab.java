package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author mahony
 * Represents an entry in the core.lab table.  
 * This is used as metadata for sequencing experiments.  
 * The name reflects the lab in which the experiment was performed (e.g. "Pugh" or "Mazzoni").
 *
 */
public class Lab implements Comparable<Lab>{
	private int dbid;
	private String name;
    
    /**
     * Creates a new <code>Lab</code> object from
     * a ResultSet that is set to a row with two
     * values:
     *  an integer that is the database id
     *  a string that is the name
     */
    public Lab(ResultSet rs) throws SQLException { 
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
		if(!(o instanceof Lab)) { return false; }
		Lab c = (Lab)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(Lab l) { return name.compareTo(l.name); }
}
