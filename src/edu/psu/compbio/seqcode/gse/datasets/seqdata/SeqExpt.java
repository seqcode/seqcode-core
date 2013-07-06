/*
 * Created as "ChipSeqExpt" on Sep 8, 2006
 */
package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

/**
 * @author tdanford 
 * @author mahony 
 * 
 * Represents a row in the seqexpts table, from the *seqdata schema.
 * 
 * 
	 create table seqexpt (
		id int(10) not null auto_increment,
		name varchar(200) not null,
		replicate varchar(200) not null,
		species int(10) not null,
		expttype int(10) not null,
		lab int(10) not null,
		exptcondition int(10) not null,
		expttarget int(10) not null,
		cellline int(10) not null,
		readtype int(10) not null,
		readlength int(5) not null,
		numreads int(12),
		collabid varchar(200),
		publicsource varchar(200),
		publicdbid varchar(200),
		fqfile varchar(500),
		exptnote longtext,
		unique (id),
		index (name, replicate, expttype, lab, exptcondition, expttarget, cellline),
		primary key (name, replicate)
	)Type=InnoDB;
 */
public class SeqExpt {
    
    private int dbid;
    private String name, replicate;
    private ExptType type;
    private Lab lab;
    private ExptTarget target;
    private ExptCondition condition;
    private CellLine cells;
    private ReadType readType;
    private Organism species;
    private int readlength, numReads;
    private String collabID, publicSource, publicDBID, fqFile, exptNote;
    
    public SeqExpt(ResultSet rs, SeqDataLoader loader) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
        replicate = rs.getString(3);
        int speciesID = rs.getInt(4);
        try {
            species = Organism.getOrganism(speciesID);
        } catch (NotFoundException e) {
            species = null;
            e.printStackTrace();
        }
        int exptTypeID = rs.getInt(5);
        int labID = rs.getInt(6);
        int condID = rs.getInt(7);
        int targetID = rs.getInt(8);
        int cellsID = rs.getInt(9);
        int readTypeID = rs.getInt(10);
        MetadataLoader mloader = loader == null ? new MetadataLoader() : loader.getMetadataLoader();
        type = mloader.loadExptType(exptTypeID);
        lab = mloader.loadLab(labID);
        condition = mloader.loadExptCondition(condID);
        target = mloader.loadExptTarget(targetID);
        cells = mloader.loadCellLine(cellsID);
        readType = mloader.loadReadType(readTypeID);

        readlength = rs.getInt(11);
        numReads = rs.getInt(12);
        collabID = rs.getString(13);
        publicSource = rs.getString(14);
        publicDBID = rs.getString(15);
        fqFile = rs.getString(16);
        exptNote = rs.getString(17);
    }
    
    public int getDBID() { return dbid; }
    public String getName() { return name; }
    public String getReplicate() { return replicate; }
    public Organism getOrganism() { return species; }
    public Lab getLab() {return lab;}
    public ExptType getExptType() {return type;}
    public ExptTarget getExptTarget() { return target; }
    public CellLine getCellLine() { return cells; }
    public ExptCondition getExptCondition() { return condition; }
    public ReadType getReadType() {return readType;}
    public int getReadLength() {return readlength;}
    public int getNumRead() {return numReads;}
    public String getCollabID() {return collabID;}
    public String getPublicSource() {return publicSource;}
    public String getPublicDBID() {return publicDBID;}
    public String getFQFile() {return fqFile;}
    public String getExptNote(){return exptNote;}
    
    public String toString() { 
    	return String.format("%s (%s)", name, replicate);
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof SeqExpt)) { return false; }
        SeqExpt c = (SeqExpt)o;
        if(dbid != c.dbid) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        return code;
    }
    
    public static PreparedStatement createLoadAll(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote from seqexpt");
    }
    
    public static PreparedStatement createLoadByDBID(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote from seqexpt where id=?");
    }
    
    public static PreparedStatement createLoadByName(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote " +
        		"from seqexpt where name=?");
    }
    
    public static PreparedStatement createLoadByLab(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote " +
        		"from seqexpt where lab=?");
    }
    
    public static PreparedStatement createLoadByCondition(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote " +
        		"from seqexpt where exptcondition=?");
    }
    
    public static PreparedStatement createLoadByTarget(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote " +
        		"from seqexpt where expttarget=?");
    }
    
    public static PreparedStatement createLoadByCellline(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote " +
        		"from seqexpt where cellline=?");
    }
    
    public static PreparedStatement createLoadByNameReplicate(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement(
        		"select id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote " +
        		"from seqexpt where name=? and replicate=?");
    }
    
    public static PreparedStatement createInsert(java.sql.Connection c) throws SQLException { 
    	String query = String.format(
                "insert into seqexpt (name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote) " +
                "values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
    	return c.prepareStatement(query);
    }
    public static PreparedStatement createInsertWithID(java.sql.Connection c) throws SQLException { 
    	String query = String.format(
                "insert into seqexpt (id, name, replicate, species, expttype, lab, exptcondition, expttarget, cellline, readtype, readlength, numreads, collabid, publicsource, publicdbid, fqfile, exptnote) " +
                "values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
    	return c.prepareStatement(query);
    }
    public static PreparedStatement createUpdateWithID(java.sql.Connection c) throws SQLException { 
    	String query = String.format(
                "update seqexpt set name=?, replicate=?, species=?, expttype=?, lab=?, exptcondition=?, expttarget=?, cellline=?, readtype=?, readlength=?, numreads=?, collabid=?, publicsource=?, publicdbid=?, fqfile=?, exptnote=? " +
                " where id=?");
    	return c.prepareStatement(query);
    }
    public static PreparedStatement createDeleteByDBID(java.sql.Connection c) throws SQLException { 
        return c.prepareStatement("delete from seqexpt where id=?");
    }
}
