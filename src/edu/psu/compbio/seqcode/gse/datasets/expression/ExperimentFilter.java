/**
 * 
 */
package edu.psu.compbio.seqcode.gse.datasets.expression;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 * @author Timothy Danford
 *
 */
public class ExperimentFilter implements Closeable {

	private MetadataLoader chipLoader;
	private ExpressionLoader loader;
	private Genome genome;
	private boolean closed;
	
	public ExperimentFilter(Genome g, ExpressionLoader l) 
		throws SQLException, UnknownRoleException {
		
		genome = g;
		loader = l;
		chipLoader = new MetadataLoader();
		closed = false;
	}
	
	public boolean isClosed() { return closed; }
	
	public void close() { 
		chipLoader.close();
		loader = null;
		closed = true;
	}
	
	public Genome getGenome() { return genome; }
	public void setGenome(Genome g) { genome = g; }
	
	public CellLine getCells(Experiment expt) 
		throws SQLException, UnknownRoleException {
		
		java.sql.Connection cxn = loader.getConnection();
		PreparedStatement ps = cxn.prepareStatement("select cells from experiment" +
				" where id=?");
		ps.setInt(1, expt.getDBID());
		ResultSet rs = ps.executeQuery();
		int cellsID = -1;
		if(rs.next()) { 
			cellsID = rs.getInt(1);
		}
		rs.close();
		ps.close();
		
		return chipLoader.loadCellLine(cellsID);
	}
	
	public Collection<Experiment> findExpts(CellLine cells, 
			ExptCondition cond) throws SQLException  {
		
        LinkedList<Integer> dbids = new LinkedList<Integer>();
        
        String tables = "experiment e";
        if(genome != null) { tables += ", probe_platform_to_genome p2g"; }
        
        String query = "select e.id from " + tables;
        		
        Vector<String> conds = new Vector<String>();
        
        if(genome != null) { 
        	int gid = genome.getDBID();
        	conds.add("e.platform=p2g.platform and p2g.genome=" + gid);
        }
        
        if(cells != null) {
        	conds.add("e.cells=" + cells.getDBID());
        }
        
        if(cond != null) { 
        	conds.add("e.condition=" + cond.getDBID());
        }
        
        for(int i = 0; i < conds.size(); i++) { 
        	if(i == 0) { 
        		query += " where ";
        	} else { 
        		query += " and ";
        	}
        	query += conds.get(i);
        }
        
        System.out.println("Final Query: \"" + query + "\"");
        
        java.sql.Connection cxn = loader.getConnection();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery(query);
        
        while(rs.next()) { 
        	int dbid = rs.getInt(1);
        	dbids.addLast(dbid);
        }
        
        rs.close();
        s.close();
        
		Collection<Experiment> expts = loader.loadExperiments(dbids);
		return expts;
	}
}
