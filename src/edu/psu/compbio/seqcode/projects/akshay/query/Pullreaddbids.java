package edu.psu.compbio.seqcode.projects.akshay.query;

import java.sql.*;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.gse.datasets.core.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseConnectionManager;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;


public class Pullreaddbids extends Pullexpttable {
	public Genome genome;
	public String queryalign = "select id from seqalignment where expt=?";
	public AlignType aligntype=null;
	/*
	 * extends to the seqalignmetn SQL table from the seqexpt SQL table
	 * genome and alignmentype are the extra options
	 */
	public Pullreaddbids(String[] command) throws SQLException,NotFoundException{
		super(command);
		MetadataLoader core = new MetadataLoader();
		try{
			
			for(int i=0; i < command.length-1; i++){
				if(command[i].matches("^--.*$") && (command[i].substring(2).matches("genome"))){
				 genome = new Genome(species,command[i+1]);
				}
				if(command[i].matches("^--.*$") && (command[i].substring(2).matches("aligntype"))){
					aligntype = core.loadAlignType(command[i+1], false, false);
				 
				}
			}
		}
		catch (NullPointerException e){
			throw new NullPointerException("Cannot use genome parameter without a specie parameter");
		}
	}
	/*
	 * Prepares and executes the SQL command and populates the resultset
	 * 
	 * Akshay: this doesn't really make sense - this method isn't called by anything. -Shaun 
	 */
	public void executeSQLAlignCommand(ResultSet rsexpt) throws SQLException, NotFoundException{
		PreparedStatement psalign=null;
		ResultSet rsalign=null;
		Connection cxn = null;
		try{
			cxn = DatabaseConnectionManager.getConnection(role);
			while(rsexpt.next()){
				boolean and = true;
				if (genome != null){queryalign += (and ? " and ": " ")+"genome = ?"; and=true;  }
				if (aligntype != null){queryalign += (and ? " and ": " ")+"aligntype = ?"; and = true;}
				psalign = cxn.prepareStatement(queryalign);
				psalign.setInt(1, rsexpt.getInt(1));
				Integer count = 2;
				if (genome != null){ psalign.setInt(count, genome.getDBID()); count += 1;}
				if (aligntype != null){ psalign.setInt(count, aligntype.getDBID() ); count += 1;}
				rsalign = psalign.executeQuery();
				getTableFromRS(rsalign);
				psalign.close();
				rsalign.close();
				queryalign="select id from seqalignment where expt=?";
			}
		}
		catch (NullPointerException e){
			e.printStackTrace();
			throw new NullPointerException("populate the rsexpt by running the executeSQLExptCommand or there are no entries in the " +
					"database for the given input parameters ");
		} finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		
	}
	
	/*
	 * iterates over the resultser and generates a table from it
	 */
	
	private void getTableFromRS(ResultSet rsalign) throws SQLException, NotFoundException{
		while(rsalign.next()){
			table.add(rsalign.getString("id"));
		}
	}
	
	/*
	 * prints the list: table
	 */
	public void printTable(){
		try{
			for(int i=0; i<table.size(); i++){
				System.out.println(table.get(i));
			}
		}
		catch (NullPointerException e){
			throw new NullPointerException("Populate the table list by calling the executeSQLAlignCommand command OR there are no entries in the " +
					"database for the given input parameters");
			
		}
		
	}

}
