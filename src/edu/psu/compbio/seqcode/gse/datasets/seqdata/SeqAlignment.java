package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import edu.psu.compbio.seqcode.gse.datasets.general.AlignType;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

import java.sql.*;
import java.util.*;

/**
 * @author tdanford
 * @author mahony
 *
 * Represents a row from the seqalignment table, in the *seqdata schemas.
 * 
 * create table seqalignment (
 * 	id int(10) not null auto_increment,
 * 	expt int(10) not null,
 * 	name varchar(200) not null,
 * 	genome int(10) not null,
 * 	permissions varchar(500),
 *  aligntype int(10) not null,
 * 	numhits int(15),
 * 	totalweight float(17,2),
 * 	aligndir varchar(200),
 * 	collabalignid varchar(200),
 * 	foreign key fk_seqalignment_expt (expt) references seqexpt(id),
 * 	index(name),
 * 	unique(id),
 * 	primary key(name)
 * );
 */
public class SeqAlignment {
	
	private int dbid;
	private SeqExpt expt;
	private String name;
	private Genome genome;
	private List<String> permissions;
	private AlignType atype;
	private int numHits;
	private double totalWeight;
	private String alignDir, collabAlignID;
	
	public SeqAlignment(ResultSet rs, SeqDataLoader loader) throws SQLException, NotFoundException { 
		dbid = rs.getInt(1);
		int exptID = rs.getInt(2);
		name = rs.getString(3);
		int genomeID = rs.getInt(4);
		int atypeID =  rs.getInt(6);
		try {
			expt = loader.loadExperiment(exptID);
			genome = Organism.findGenome(genomeID);
			MetadataLoader mloader = loader == null ? new MetadataLoader() : loader.getMetadataLoader();
			atype = mloader.loadAlignType(atypeID);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		permissions = Arrays.asList(rs.getString(5).split(";"));
		numHits = rs.getInt(7);
		totalWeight = rs.getFloat(8);
		alignDir = rs.getString(9);
		collabAlignID = rs.getString(10);
	}
	
	public SeqAlignment(ResultSet rs, SeqExpt ex, AlignType a) throws SQLException { 
		dbid = rs.getInt(1);
		expt = ex;  // exptID is still in position #2, we just don't need it.
		name = rs.getString(3);
		int genomeID = rs.getInt(4);
		atype =  a;
		try {
			genome = Organism.findGenome(genomeID);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		permissions = Arrays.asList(rs.getString(5).split(";"));
		numHits = rs.getInt(7);
		totalWeight = rs.getFloat(8);
		alignDir = rs.getString(9);
		collabAlignID = rs.getString(10);
	}
	
	public int getDBID() { return dbid; }
	public String getName() { return name; }
	public SeqExpt getExpt() { return expt; }
	public AlignType getAlignType() {return atype;}
	public Genome getGenome() { return genome; }
	public List<String> getPermissions() { return permissions; }
	public int getNumHits() {return numHits;}
	public double getTotalWeight() {return totalWeight;}
	public String getAlignDir() {return alignDir;}
	public String getCollabAlignID() {return collabAlignID;}
	
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public String toString() { 
		return String.format("\"%s\" %s -> %s", name, expt.getName(), genome.getVersion());
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof SeqAlignment)) { return false; }
		SeqAlignment a = (SeqAlignment)o;
		if(dbid != a.dbid) { return false; }
		return true;
	}

	public static PreparedStatement createLoadByIDStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid from seqalignment where id=?";
		return c.prepareStatement(query);
	}
	
	public static PreparedStatement createLoadByNameAndExptStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid from seqalignment where name=? and expt=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createLoadAllByExptStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid from seqalignment where expt=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createLoadAllByGenomeStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid from seqalignment where genome=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createLoadAllByExptAndGenomeStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid from seqalignment where expt=? and genome=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createInsertStatement(java.sql.Connection c) throws SQLException { 
		String query = String.format(
				"insert into seqalignment (expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid) " +
				"values (?, ?, ?, ?, ?, ?, ?, ?, ?)");
		return c.prepareStatement(query);
	}
	
	public static PreparedStatement createInsertStatementWithID(java.sql.Connection c) throws SQLException { 
		String query = String.format(
				"insert into seqalignment (id, expt, name, genome, permissions, aligntype, numhits, totalweight, aligndir, collabalignid) " +
				"values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
		return c.prepareStatement(query);
	}
	
	public static PreparedStatement updateHitsAndWeights(java.sql.Connection c) throws SQLException { 
		String query = "update seqalignment set numhits=?, totalweight=? where id=?";
		return c.prepareStatement(query);
	}
	
}
