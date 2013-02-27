package edu.psu.compbio.seqcode.gse.datasets.chipchip;

import java.sql.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.locators.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class ChipChipFilter implements Closeable {
	
	private Genome genome;
	private java.sql.Connection cxn;
	
	public ChipChipFilter(Genome g) throws SQLException, UnknownRoleException {
		genome = g;
		cxn = DatabaseFactory.getConnection(ExptLocator.dbRole);
	}
	
	public Pair<ExptTarget,ExptTarget> getFactorData(MetadataLoader loader, ChipChipLocator loc) throws SQLException { 
		Pair<ExptTarget,ExptTarget> p = null;
		String query = "select e.factorone, e.factortwo from experiment e where e.name=? and e.version=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, loc.getNameVersion().name);
		ps.setString(2, loc.getNameVersion().version);

		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			int c1 = rs.getInt(1), c2 = rs.getInt(2);
			ExptTarget cond1 = loader.loadExptTarget(c1);
			ExptTarget cond2 = loader.loadExptTarget(c2);
			p = new Pair<ExptTarget,ExptTarget>(cond1, cond2);
		} 
		rs.close();
		ps.close();

		return p;		
	}
	
	public Pair<ExptCondition,ExptCondition> getConditionData(MetadataLoader loader, ChipChipLocator loc) throws SQLException { 
		Pair<ExptCondition,ExptCondition> p = null;
		String query = "select e.conditionone, e.conditiontwo from experiment e where e.name=? and e.version=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, loc.getNameVersion().name);
		ps.setString(2, loc.getNameVersion().version);

		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			int c1 = rs.getInt(1), c2 = rs.getInt(2);
			ExptCondition cond1 = loader.loadExptCondition(c1);
			ExptCondition cond2 = loader.loadExptCondition(c2);
			p = new Pair<ExptCondition,ExptCondition>(cond1, cond2);
		} else { 
		    System.err.println("Couldn't find any values for locator: " + loc.toString());
        }
		rs.close();
		ps.close();

		return p;		
	}
	
	public Pair<CellLine,CellLine> getCellsData(MetadataLoader loader, ChipChipLocator loc) throws SQLException {
		Pair<CellLine,CellLine> p = null;
		String query = "select e.cellsone, e.cellstwo from experiment e where e.name=? and e.version=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, loc.getNameVersion().name);
		ps.setString(2, loc.getNameVersion().version);

		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			int c1 = rs.getInt(1), c2 = rs.getInt(2);
			CellLine cells1 = loader.loadCellLine(c1);
			CellLine cells2 = loader.loadCellLine(c2);
			p = new Pair<CellLine,CellLine>(cells1, cells2);
		} 
		rs.close();
		ps.close();

		return p;
	}

	public Collection<ExptLocator> 
		findBinding(CellLine cells, ExptCondition cond, ExptTarget factor) throws SQLException {
		
		LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>();
		Statement s = cxn.createStatement();
		
		String agilentQuery = "select e.name, e.version from experiment e, " +
				"exptToGenome e2g where e.id=e2g.experiment and e2g.genome=" + 
				genome.getDBID();
		
		if(cells != null) { 
			int dbid = cells.getDBID();
			agilentQuery += " and (e.cellsone=" + dbid + " or e.cellstwo=" + dbid + ")";
		}
		
		if(cond != null) { 
			int dbid = cond.getDBID();
			agilentQuery += " and (e.conditionone=" + dbid + " or e.conditiontwo=" + dbid + ")";
		}

		if(factor != null) { 
			int dbid = factor.getDBID();
			agilentQuery += " and (e.factorone=" + dbid + " or e.factortwo=" + dbid + ")";
		}
		
		ResultSet rs = s.executeQuery(agilentQuery);
		while(rs.next()) { 
			String n = rs.getString(1), v = rs.getString(2);
			ChipChipLocator loc = new ChipChipLocator(genome, n, v);
			locs.addLast(loc);
		}
		rs.close();
		
		s.close();
		return locs;
	}

	public void close() {
		DatabaseFactory.freeConnection(cxn);
		cxn=null;
	}

	public boolean isClosed() {
		return cxn==null;
	}
}
