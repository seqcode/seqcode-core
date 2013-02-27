package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author mahony
 * Represents an entry in the core.expttype table.  
 * This is used as metadata for sequencing experiments.  
 * The name reflects the type of sequencing experiment (e.g. "CHIPSEQ", "RNASEQ", "CHIPEXO", etc)
 *
 */
public class ExptType implements Comparable<ExptType>{
	private int dbid;
	private String name;
    
    /**
     * Creates a new <code>ExptType</code> object from
     * a ResultSet that is set to a row with two
     * values:
     *  an integer that is the database id
     *  a string that is the name
     */
    public ExptType(ResultSet rs) throws SQLException { 
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
		if(!(o instanceof ExptType)) { return false; }
		ExptType c = (ExptType)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(ExptType e) { return name.compareTo(e.name); }
}
