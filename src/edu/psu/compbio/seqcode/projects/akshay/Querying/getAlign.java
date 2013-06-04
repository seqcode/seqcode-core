package edu.psu.compbio.seqcode.projects.akshay.Querying;
import java.sql.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;

public class getAlign extends getExpt{
	public Genome genome;
	public String queryalign = "select id from seqalignment where expt=?";
	public AlignType aligntype=null;
	public PreparedStatement psalign=null;
	public ResultSet rsalign=null;
	
	public getAlign(String[] command) throws SQLException,NotFoundException{
		super(command);
		try{
			for(int i=0; i < command.length-1; i++){
				if(command[i].matches("^--.*$") && (command[i].substring(2).matches("genome"))){
				 genome = new Genome(species.getName(),command[i+1]);
				}
				if(command[i].matches("^--.*$") && (command[i].substring(2).matches("aligntype"))){
					aligntype = core.getAlignType(command[i+1]);
				 
				}
			}
		}
		catch (NullPointerException e){
			throw new NullPointerException("Cannot use genome parameter without a specie parameter");
		}
	}
	
	public void executeSQLAlignCommand() throws SQLException, NotFoundException{
		try{
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
				getTableFromRS();
				psalign=null;
				rsalign=null;
				queryalign="select id from seqalignment where expt=?";
				
			}
		}
		catch (NullPointerException e){
			e.printStackTrace();
			throw new NullPointerException("populate the rsexpt by running the executeSQLExptCommand or there are no entries in the " +
					"database for the given input parameters ");
		}
		
	}
	
	public void getTableFromRS() throws SQLException, NotFoundException{
		while(rsalign.next()){
			table.add(rsalign.getString("id"));
		}
	}
	
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
	
	public void closeConnection() throws SQLException{
		core.close();
		if (cxn != null){
			DatabaseFactory.freeConnection(cxn);
	        cxn = null;
		}
	if (psexpt != null){psexpt.close();}
	if (rsexpt != null){rsexpt.close();}
	if (psalign != null){psalign.close();}
	if (rsalign != null){rsalign.close();}
	
	
	}
	
	
		
}
		

