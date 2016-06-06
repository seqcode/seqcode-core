/**
 * 
 */
package org.seqcode.data.core;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author tdanford
 * Represents an entry in the core.expttarget table.  
 * This is used as metadata for sequencing experiments.  
 * The name might describe describe the target of the antibody in ChIP-seq experiments, or the type of RNA targeted by RNA-seq experiments. 
 */
public class ExptTarget implements Comparable<ExptTarget> {
    
	private int dbid;
	private String name;
    
    public ExptTarget(ResultSet rs) throws SQLException { 
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
		if(!(o instanceof ExptTarget)) { return false; }
		ExptTarget c = (ExptTarget)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
    
    public int compareTo(ExptTarget f) { return name.compareTo(f.name); }

    
    public static PreparedStatement prepareLoadExperimentsByGenome(java.sql.Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select e.id, e.name, e.version from experiment e, exptToGenome e2g " +
                "where e.id=e2g.experiment and e2g.genome=? and e.active=1 and (e.factorone=? or e.factortwo=?)");
    }    
}
