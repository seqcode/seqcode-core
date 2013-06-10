/**
 * 
 */
package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collection;
import java.util.LinkedList;

import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

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

    public Collection<ExptLocator> loadLocators(java.sql.Connection cxn, Genome g) throws SQLException {
        LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>();
        PreparedStatement ps = prepareLoadExperimentsByGenome(cxn);
        
        ps.setInt(1, g.getDBID());
        ps.setInt(2, dbid); ps.setInt(3, dbid);
        
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            String name = rs.getString(2);
            String version = rs.getString(3);
            ChipChipLocator loc = new ChipChipLocator(g, name, version);
            locs.addLast(loc);
        }
        
        rs.close();
        ps.close();
        return locs;
    }
    
    public static PreparedStatement prepareLoadExperimentsByGenome(java.sql.Connection cxn) throws SQLException { 
        return cxn.prepareStatement("select e.id, e.name, e.version from experiment e, exptToGenome e2g " +
                "where e.id=e2g.experiment and e2g.genome=? and e.active=1 and (e.factorone=? or e.factortwo=?)");
    }    
}
