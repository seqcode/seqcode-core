package edu.psu.compbio.seqcode.projects.akshay.query;

import java.sql.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;

public class Pullexpttable {
	
	public ExptType expttype=null;
	public Lab lab=null;
	public ExptCondition exptcondition=null;
	public ExptTarget expttarget=null;
	public CellLine celline=null;
	public ReadType readtype=null;
	public Organism species=null;
	public List<String> table= new ArrayList<String>();
	
	public java.sql.Connection cxn;
	public static final String role = "seqdata";
	public String queryexpt = "select id, name, replicate, species, expttype," +
			" exptcondition, lab, expttarget, cellline, readtype from seqexpt";
	public ResultSet rsexpt=null;
	
	/*public getExpt(){}*/
	public MetadataLoader core = new MetadataLoader();
	public PreparedStatement psexpt=null;
	
	
	public Pullexpttable(String[] command) throws SQLException{
		
		try{
		
		cxn = DatabaseFactory.getConnection(role);
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
	
	public void executeSQLExptCommand() throws SQLException{
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
		System.out.println(queryexpt);
		System.out.println(expttarget.getDBID());
		System.out.println(celline.getDBID());
		Integer count = 1;
		
		if (lab != null){ psexpt.setInt(count, lab.getDBID()); count += 1;}
		if (expttype != null){ psexpt.setInt(count, expttype.getDBID()); count += 1;}
		if (expttarget != null){psexpt.setInt(count, expttarget.getDBID()); count += 1;}
		if (exptcondition != null) {psexpt.setInt(count, exptcondition.getDBID()); count += 1;}
		if (celline != null){psexpt.setInt(count, celline.getDBID()); count += 1;}
		if (readtype != null){psexpt.setInt(count, readtype.getDBID()); count+=1;}
		if (species != null){ psexpt.setInt(count, species.getDBID()); count +=1;}
	
		
		rsexpt = psexpt.executeQuery();
		/*psexpt.clearParameters();*/
		
			
		}
	
	public void getTableFromRS() throws SQLException, NotFoundException{ 
		try{
			MetadataLoader trackback = new MetadataLoader();
			Organism org=null;
		
			while(rsexpt.next()){
				System.out.println(rsexpt.getString(1));
				org = new Organism(Integer.parseInt(rsexpt.getString("species")));
				
				table.add(rsexpt.getString(1)+"\t"+rsexpt.getString(1) + "\t"+rsexpt.getString(2)+"\t"+org.getName()+"\t"+
						trackback.loadExptType(Integer.parseInt(rsexpt.getString("expttype"))).getName()+"\t"+trackback.loadLab(Integer.parseInt(rsexpt.getString("lab"))).getName()+"\t"+
						trackback.loadExptCondition(Integer.parseInt(rsexpt.getString("exptcondition"))).getName()+"\t"+trackback.loadExptTarget(Integer.parseInt(rsexpt.getString("expttarget"))).getName()+"\t"+
						trackback.loadCellLine(Integer.parseInt(rsexpt.getString("cellline"))).getName()+"\t"+trackback.loadReadType(Integer.parseInt(rsexpt.getString("readtype"))).getName());
			}
			trackback.close();
		}
		catch (NullPointerException e){
			throw new NullPointerException("Populate the ResultSet object by calling the executeSQLCommand method");
		}
		
	}
	
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
	
	public void closeConnection() throws SQLException{
		core.close();
		if (cxn != null){
		DatabaseFactory.freeConnection(cxn);
        cxn = null;
        }
		if (psexpt != null){psexpt.close();}
		if (rsexpt != null){rsexpt.close();}
		
	}

}
