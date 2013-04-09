package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author mahony
 * Represents an entry in the core.aligntype table.  
 * This is used as metadata for sequencing experiments.  
 * The name reflects the type of sequencing alignment, which is "SINGLE" or "PAIRED" for now.
 * Looks redundant with readtype - but you can have a "PAIRED" read that is "SINGLE" aligned,
 * and there could be more variation in the future.
 *
 */
public class AlignType implements Comparable<AlignType>{
	private int dbid;
	private String name;
    
    /**
     * Creates a new <code>AlignType</code> object from
     * a ResultSet that is set to a row with two
     * values:
     *  an integer that is the database id
     *  a string that is the name
     */
    public AlignType(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
    }
	
	public String getName() { return name; }
	public int getDBID() { return dbid; }
    
    public String toString() { return name; }
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof AlignType)) { return false; }
		AlignType c = (AlignType)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(AlignType r) { return name.compareTo(r.name); }
}
