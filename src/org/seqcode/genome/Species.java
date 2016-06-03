package org.seqcode.genome;

import java.util.*;


import java.sql.*;

import org.seqcode.gse.utils.*;
import org.seqcode.gse.utils.database.DatabaseConnectionManager;
import org.seqcode.gse.utils.database.DatabaseException;
import org.seqcode.gse.utils.database.Sequence;
import org.seqcode.gse.utils.database.UnknownRoleException;

/**
 * Species represents a species in our database
 */
public class Species{

    //static cache of all Species
    private static Map<String,Species> organisms = new HashMap<String,Species>();
    private static Map<Integer, Species> organismsids = new HashMap<Integer,Species>();

    //Variables that define this Species
    private String species;
    private int dbid;
    
    /**
     * Constructor: no db lookup needed
     * @param id : use -1 for fake organisms (i.e. not in db)
     * @param species
     * @throws NotFoundException
     */
    public Species(int id, String species){
    	this.dbid = id;
    	this.species = species;
    	
    }
    
    /**
     * Constructor: look up species name in core database
     * @param species
     * @throws NotFoundException
     */
    public Species(String species) throws NotFoundException {
    	this.species = species;
    	if(organisms.isEmpty()){
    		getAllSpecies(false);
        }
        if (organisms.containsKey(species)) {
            this.dbid = organisms.get(species).getDBID();
        }else{ //Try looking up in the database, in case it was put there since organismids was populated
	        Connection cxn = null;
	        Statement stmt = null;
	        ResultSet rs = null;
	        try {
	            cxn = DatabaseConnectionManager.getConnection("core");
	            stmt = cxn.createStatement();
	            rs = stmt.executeQuery("select id from species where name = '" + species + "'");
	            
	            if (rs.next()) {
	                this.dbid = rs.getInt(1);
	            } else {
	                throw new NotFoundException("Couldn't find " + species);
	            } 
	
	        } catch (SQLException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("Couldn't find " + species + ": " + ex.toString(), ex);
	        } catch (UnknownRoleException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("Couldn't connect with role core", ex);
	        } finally {
	        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
		        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
	        	if (cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	        }
        }
    }


    /**
     * Constructor: look up dbid in core database
     * @param speciesID
     * @throws NotFoundException
     */
    public Species(int speciesID) throws NotFoundException {
        this.dbid = speciesID;
        if(organisms.isEmpty()){
        	getAllSpecies(false);
        }
        if (organisms.containsKey(dbid)) {
            this.species = organisms.get(dbid).getName();
            
        }else{ //Try looking up in the database, in case it was put there since organismids was populated
	        java.sql.Connection cxn = null;
	        Statement stmt = null;
	        ResultSet rs = null;
	        try {
	            cxn = DatabaseConnectionManager.getConnection("core");
	            stmt = cxn.createStatement();
	            rs = stmt.executeQuery("select name from species where id=" + dbid);
	            if (rs.next()) {
	                species = rs.getString(1);
	            } else {
	                throw new NotFoundException("Couldn't find " + dbid);
	            }
	        } catch (SQLException ex) {
	            throw new DatabaseException("Couldn't find " + dbid + ": " + ex.toString(), ex);
	        } catch (UnknownRoleException ex) {
	            throw new DatabaseException("Couldn't connect with role core", ex);
	        } finally {
	        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
		        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
	        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	        }
        }
        
    }



    /**
     * Accessor for species name
     */
    public String getName() {
        return species;
    }


    /**
     * Accessor for database ID
     */
    public int getDBID() {
        return dbid;
    }

	/**
	 * Returns a list of genomes for this species
	 * @return
	 */
    public Collection<Genome> getGenomes(){
    	try {
			return Genome.getAllGenomesBySpecies(this);
		} catch (NotFoundException e) {
			e.printStackTrace();			
		}
    	return null;
    }
    
    /**
	 * Returns a list of genome names for this species
	 * @return
	 */
    public Collection<String> getGenomeNames(){
    	List<String> names = new ArrayList<String>();
    	for(Genome g : getGenomes()){
    		names.add(g.getVersion()); 
    	}
    	return names;
    }
    
	/**
	 * Load all Organisms from database 
	 * @return
	 */
    public static Collection<Species> getAllSpecies(boolean forceRefreshFromDB){
    	List<Species> orgs = new ArrayList<Species>();
    	
    	if(organisms.isEmpty() || forceRefreshFromDB){
    		organisms.clear(); organismsids.clear();
	    	Connection cxn = null;
	        Statement stmt = null;
	        ResultSet rs = null;
	        try {
	            cxn = DatabaseConnectionManager.getConnection("core");
	            stmt = cxn.createStatement();
	            rs = stmt.executeQuery("select id, name from species");
	            while(rs.next()) { 
	            	Species org = new Species(rs.getInt(1), rs.getString(2));
	            	orgs.add(org);
	            	organisms.put(org.getName(), org);
	                organismsids.put(org.getDBID(), org);
	            }
	
	        } catch (SQLException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("mySQL error: " + ex.toString(), ex);
	        } catch (UnknownRoleException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("Couldn't connect with role core", ex);
	        } finally {
	        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
		        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
	        	if (cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	        }
    	}else{
    		orgs.addAll(organisms.values());
    	}
	    return orgs;
    }

    /**
     * Get a specific Species by name
     * @param species
     * @return
     * @throws NotFoundException
     */
    public static Species getSpecies(String species) throws NotFoundException {
        if(organisms.isEmpty()){
        	getAllSpecies(false);
        }
    	if (organisms.containsKey(species)) {
            return organisms.get(species);
        } else {
            Species o = new Species(species);
            organisms.put(species, o);
            organismsids.put(o.getDBID(), o);
            return o;
        }
    }


    /**
     * Get a specific Species by ID
     * @param species
     * @return
     * @throws NotFoundException
     */
    public static Species getSpecies(int species) throws NotFoundException {
    	if(organisms.isEmpty()){
        	getAllSpecies(false);
        }
        if (organismsids.containsKey(species)) {
            return organismsids.get(species);
        } else {
            Species o = new Species(species);
            organisms.put(o.getName(), o);
            organismsids.put(o.getDBID(), o);
            return o;
        }
    }
    
    /**
     * Get all Species names from the database
     * @return
     */
    public static Collection<String> getAllSpeciesNames(boolean forceRefreshFromDB) {
    	if(organisms.isEmpty() || forceRefreshFromDB){
        	getAllSpecies(forceRefreshFromDB);
        }
    	return organisms.keySet();
    }

    /**
     * Insert a new Species into the database
     * @param species
     * @throws SQLException
     */
    public static void insertSpecies(String species) throws SQLException {
        java.sql.Connection cxn = null;
        Statement stmt = null;
        
        try {
            cxn = DatabaseConnectionManager.getConnection("core");
            stmt = cxn.createStatement();
            String nextIdString = Sequence.getInsertSQL(cxn, "species_id");
            String sql = String.format("insert into species values (%s, '%s')", nextIdString, species);
            System.out.println(sql);
            stmt.executeUpdate(sql);
            
            Species o = new Species(new Integer(nextIdString), species);
            organisms.put(species, o);
            organismsids.put(o.getDBID(), o);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role core", ex);
        } catch (SQLException se) {
            throw se;
        } finally {
	        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
        }
    }


}
