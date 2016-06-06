package org.seqcode.data.core;

import java.sql.*;


/**
 * @author tdanford
 *
 * Represents an entry in the core.cellline table.  This is used as metadata
 * for sequencing experiments.  The name might be something
 * general (and not very useful) like 'S288C' or a cell line name or strain number.
 */
public class CellLine implements Comparable<CellLine> {
	
	private int dbid;
	private String name;
    
    /**
     * Creates a new <code>CellLine</code> object from
     * a ResultSet that is set to a row with two
     * values:
     *  an integer that is the database id
     *  a string that is the name
     */
    public CellLine(ResultSet rs) throws SQLException { 
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
		if(!(o instanceof CellLine)) { return false; }
		CellLine c = (CellLine)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(CellLine c) { return name.compareTo(c.name); }
    
    
    public static PreparedStatement prepareLoadExperimentsByGenome(java.sql.Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select e.id, e.name, e.version from experiment e, exptToGenome e2g " +
                "where e.id=e2g.experiment and e2g.genome=? and e.active=1 and (e.cellsone=? or e.cellstwo=?)");
    }
    
    
}
