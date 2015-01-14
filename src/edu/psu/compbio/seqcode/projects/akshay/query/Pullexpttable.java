package edu.psu.compbio.seqcode.projects.akshay.query;

import java.sql.*;
import java.util.*;

import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseConnectionManager;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;
import edu.psu.compbio.seqcode.gse.datasets.core.*;

public class Pullexpttable {
	
	public ExptType expttype=null;
	public Lab lab=null;
	public ExptCondition exptcondition=null;
	public ExptTarget expttarget=null;
	public CellLine celline=null;
	public ReadType readtype=null;
	public Organism species=null;
	public List<String> table= new ArrayList<String>();
	
	public static final String role = "seqdata";
	public String queryexpt = "select id, name, replicate, species, expttype," +
			" exptcondition, lab, expttarget, cellline, readtype from seqexpt";
	
	/* * command is a list of strings that is generally returned by a scanner class
	 * The following constructor class connects to the core SQL data base and initialises the different classes (celline, expttarget .. etc)
	 */
	
	public Pullexpttable(String[] command) throws SQLException{
		MetadataLoader core = new MetadataLoader();
		try{
				
			 for(int i=0; i < command.length-1; i++){
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("expttype"))){
					 expttype = core.getExptType(command[i+1]);
					
				 }
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("lab"))){
					 lab = core.getLab(command[i+1]);
					 
				 }
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("exptcondition"))){
					 exptcondition = core.getExptCondition(command[i+1]);
					 
				 }
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("expttarget"))){
					 expttarget = core.getExptTarget(command[i+1]);
					 
				 }
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("cellline"))){
					 celline = core.getCellLine(command[i+1]);
					 
				 }
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("readtype"))){
					 readtype = core.getReadType(command[i+1]);
				 }
				 if(command[i].matches("^--.*$") && (command[i].substring(2).matches("species"))){
					 species = new Organism(command[i+1]);
				 }
			}
		}
		catch(NotFoundException e){
			e.printStackTrace();
		}
		catch(UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }
	}
	
	/* * Prepares and executes the SQL command based on the input options
	 * Also populates the resultset object
	 */
	
	public void executeSQLExptCommand() throws SQLException{
		Connection cxn = null;
		PreparedStatement psexpt=null;
		MetadataLoader core = new MetadataLoader();
		try{
			cxn = DatabaseConnectionManager.getConnection(role);
			
			if((lab != null) || (expttype != null) || (expttarget != null) || (exptcondition != null) || (celline != null) || (readtype != null) || (species != null)){
				queryexpt += " where";
			}
			boolean and = false;
			if (lab != null){queryexpt += " lab = ?"; and=true;  }
			if (expttype != null){queryexpt += (and ? " and ": " ")+"expttype = ?"; and = true;}
			if (expttarget != null){queryexpt += (and ? " and ": " ")+"expttarget = ?"; and =true;}
			if (exptcondition != null){queryexpt += (and ? " and ": " ")+"exptcondition = ?"; and = true;}
			if (celline != null){queryexpt += (and ? " and ": " ")+"cellline = ?"; and = true;}
			if (readtype != null){queryexpt += (and ? " and ": " ")+"readtype = ?"; and = true;}
			if (species != null){queryexpt += (and ? " and ": " ")+"species = ?"; and = true;}
			
			psexpt = cxn.prepareStatement(queryexpt);
			
			Integer count = 1;
			
			if (lab != null){ psexpt.setInt(count, lab.getDBID()); count += 1;}
			if (expttype != null){ psexpt.setInt(count, expttype.getDBID()); count += 1;}
			if (expttarget != null){psexpt.setInt(count, expttarget.getDBID()); count += 1;}
			if (exptcondition != null) {psexpt.setInt(count, exptcondition.getDBID()); count += 1;}
			if (celline != null){psexpt.setInt(count, celline.getDBID()); count += 1;}
			if (readtype != null){psexpt.setInt(count, readtype.getDBID()); count+=1;}
			if (species != null){ psexpt.setInt(count, species.getDBID()); count +=1;}
		
			ResultSet rsexpt = psexpt.executeQuery();
			getTableFromRS(rsexpt, core);
			rsexpt.close();
			psexpt.close();
			
		} catch (NotFoundException e) {
			e.printStackTrace();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
	}
	
	/* * Iterates the resultset object. 
	 *  populates the List: table
	 */
	private void getTableFromRS(ResultSet rsexpt, MetadataLoader trackback) throws SQLException, NotFoundException{ 
		try{
			Organism org=null;
		
			while(rsexpt.next()){
				System.out.println(rsexpt.getString(1));
				org = new Organism(Integer.parseInt(rsexpt.getString("species")));
				
				table.add(rsexpt.getString(1)+ "\t"+rsexpt.getString(2)+"\t"+ rsexpt.getString("replicate")+"\t"+org.getName()+"\t"+
						trackback.loadExptType(Integer.parseInt(rsexpt.getString("expttype"))).getName()+"\t"+trackback.loadLab(Integer.parseInt(rsexpt.getString("lab"))).getName()+"\t"+
						trackback.loadExptCondition(Integer.parseInt(rsexpt.getString("exptcondition"))).getName()+"\t"+trackback.loadExptTarget(Integer.parseInt(rsexpt.getString("expttarget"))).getName()+"\t"+
						trackback.loadCellLine(Integer.parseInt(rsexpt.getString("cellline"))).getName()+"\t"+trackback.loadReadType(Integer.parseInt(rsexpt.getString("readtype"))).getName());
			}
		}
		catch (NullPointerException e){
			throw new NullPointerException("Populate the ResultSet object by calling the executeSQLCommand method");
		}
		
	}
	
	
	/* * Prints the table
	 *  Make sure to populate the table before calling this method
	 */
	public void printTable(){
		try{
			for(int i=0; i<table.size(); i++){
				System.out.println(table.get(i));
			}
		}
		catch (NullPointerException e){
			throw new NullPointerException("Populate the table list by calling the getTableFromRS command OR there are no entries in the " +
					"database for the given input parameters");
			
		}
	}

}
