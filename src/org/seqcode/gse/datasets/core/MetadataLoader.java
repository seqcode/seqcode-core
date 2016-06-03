package org.seqcode.gse.datasets.core;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;

import org.seqcode.gse.utils.database.DatabaseConnectionManager;
import org.seqcode.gse.utils.database.DatabaseException;
import org.seqcode.gse.utils.database.UnknownRoleException;


/**
 * MetadataLoader is the interface for interacting with metadata entries in the core database.
 * 
 * Methods in this class make use of a simple caching approach. Use the "forceDatabaseRefresh" flag in each method to 
 * ignore the cache and to re-query the database entries (useful if you expect that the database may have been updated
 * between queries). 
 * 
 * @author mahony
 */

public class MetadataLoader{
    
    public static final String role = "core";

    private boolean allLabsLoaded=false, allCellsLoaded=false, allCondsLoaded=false, allTargetsLoaded=false, 
    		allExptTypesLoaded=false, allReadTypesLoaded=false, allAlignTypesLoaded=false, allSeqDataUsersLoaded=false;
    
    private Map<String,Lab> labNames;
    private Map<String,CellLine> cellNames;
    private Map<String,ExptCondition> condNames;
    private Map<String,ExptTarget> targetNames;
    private Map<String,ExptType> exptTypeNames;
    private Map<String,ReadType> readTypeNames;
    private Map<String,AlignType> alignTypeNames;
    private Map<String,SeqDataUser> seqDataUserNames;
	
    private Map<Integer,Lab> labIDs;
    private Map<Integer,CellLine> cellIDs;
    private Map<Integer,ExptCondition> condIDs;
    private Map<Integer,ExptTarget> targetIDs;
    private Map<Integer,ExptType> exptTypeIDs;
    private Map<Integer,ReadType> readTypeIDs;
    private Map<Integer,AlignType> alignTypeIDs;
    private Map<Integer,SeqDataUser> seqDataUserIDs;
    	
    public MetadataLoader() throws SQLException {this(false);}
    public MetadataLoader(boolean cacheAll) throws SQLException { 
        
        labNames = new HashMap<String,Lab>();
        labIDs = new HashMap<Integer,Lab>();
        cellNames = new HashMap<String,CellLine>();
        cellIDs = new HashMap<Integer,CellLine>();
        condNames = new HashMap<String,ExptCondition>();
        condIDs = new HashMap<Integer,ExptCondition>();
        targetNames = new HashMap<String,ExptTarget>();
        targetIDs = new HashMap<Integer,ExptTarget>();
        exptTypeNames = new HashMap<String,ExptType>();
        exptTypeIDs = new HashMap<Integer,ExptType>();
        readTypeNames = new HashMap<String,ReadType>();
        readTypeIDs = new HashMap<Integer,ReadType>();
        alignTypeNames = new HashMap<String,AlignType>();
        alignTypeIDs = new HashMap<Integer,AlignType>();
        seqDataUserNames = new HashMap<String,SeqDataUser>();
        seqDataUserIDs = new HashMap<Integer,SeqDataUser>();
     
        if(cacheAll)
        	cacheAllMetadata();
    }
	
    /**
     * Load all metadata tables into the cache Maps
     * (forces update of all tables if they already exist).
     * Performs all queries with one connection establishment for efficiency.
     */
    public void cacheAllMetadata() throws SQLException {
    	Connection cxn = null;
        PreparedStatement ps=null;
        ResultSet rs = null;
		try {
			//Load all Labs
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from lab");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	        	Lab l = new Lab(rs);
	            labNames.put(l.getName(), l);
	            labIDs.put(l.getDBID(),l);
	        }
	        
	        //Load all Celllines
	        ps = cxn.prepareStatement("select id, name from cellline");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	            CellLine c = new CellLine(rs);
	            cellNames.put(c.getName(), c);
	            cellIDs.put(c.getDBID(),c);
	        }
	        
	        //Load all ExptConditions
	        ps = cxn.prepareStatement("select id, name from exptcondition");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	            ExptCondition c = new ExptCondition(rs);
	            condNames.put(c.getName(), c);
	            condIDs.put(c.getDBID(),c);
	        }
	
	        //Load all ExptTargets
	        ps = cxn.prepareStatement("select id, name from expttarget");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	            ExptTarget f = new ExptTarget(rs);
	            targetNames.put(f.getName(), f);
	            targetIDs.put(f.getDBID(),f);
	        }
	        
	        //Load all ExptTypes
	        ps = cxn.prepareStatement("select id, name from expttype");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	        	ExptType e = new ExptType(rs);
	            exptTypeNames.put(e.getName(), e);
	            exptTypeIDs.put(e.getDBID(),e);
	        }
	        
	        //Load all ReadTypes
	        ps = cxn.prepareStatement("select id, name from readtype");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	        	ReadType r = new ReadType(rs);
	            readTypeNames.put(r.getName(), r);
	            readTypeIDs.put(r.getDBID(),r);
	        }
	        
	        //Load all AlignTypes
	        ps = cxn.prepareStatement("select id, name from aligntype");
	        rs = ps.executeQuery();
	        while(rs.next()) { 
	        	AlignType a = new AlignType(rs);
	            alignTypeNames.put(a.getName(), a);
	            alignTypeIDs.put(a.getDBID(),a);
	        }
	        
	        //Load all SeqDataUsers
	        ps = cxn.prepareStatement("select id, name, admin from seqdatauser");
			rs = ps.executeQuery();
			while(rs.next()) { 
				SeqDataUser a = new SeqDataUser(rs);
				seqDataUserNames.put(a.getName(), a);
				seqDataUserIDs.put(a.getDBID(),a);
			}
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allLabsLoaded = true;
		allCellsLoaded = true;
		allCondsLoaded = true;
    	allTargetsLoaded = true;
    	allExptTypesLoaded = true;
    	allReadTypesLoaded = true;
    	allAlignTypesLoaded = true;
    	allSeqDataUsersLoaded = true;
    }
    
    
    //////////////////
    // Lab stuff
    //////////////////
    
    /**
     * Load single Lab by name
     * 
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Lab loadLab(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(labNames.containsKey(name) && !forceDatabaseRefresh) { return labNames.get(name); }
    	
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		Lab l = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from lab where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
                l = new Lab(rs);
                
                if(!labIDs.containsKey(l.getDBID())) { 
                    labIDs.put(l.getDBID(), l);
                    labNames.put(l.getName(), l);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(l == null && insertIfNone){
			int id = insertLab(name);
			return loadLab(id, forceDatabaseRefresh);
		}else{
			return l;
		}
    }
    
    /**
     * Load single Lab by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Lab loadLab(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(labIDs.containsKey(dbid) && !forceDatabaseRefresh) { return labIDs.get(dbid); }

        Lab l = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from lab where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            if(rs.next()) { 
                l = new Lab(rs);
                labIDs.put(dbid, l);
                labNames.put(l.getName(), l);
            } else {
                throw new IllegalArgumentException("Unknown Lab DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return l;
    }
    
    /**
     * Load a collection of Labs by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<Lab> loadLabs(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<Lab> values = new LinkedList<Lab>();
        for(int dbid : dbids) { values.addLast(loadLab(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all Labs
     * 
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<Lab> loadAllLabs(boolean forceDatabaseRefresh) throws SQLException {
        if(allLabsLoaded && !forceDatabaseRefresh){ return labNames.values();}
        if(forceDatabaseRefresh){labNames.clear(); labIDs.clear();}
        
    	HashSet<Lab> values = new HashSet<Lab>();
        Connection cxn = null;
        PreparedStatement loadAllLabs=null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            loadAllLabs = cxn.prepareStatement("select id, name from lab");
	        rs = loadAllLabs.executeQuery();
	
	        while(rs.next()) { 
	        	Lab l = new Lab(rs);
	            values.add(l);
	            labNames.put(l.getName(), l);
	            labIDs.put(l.getDBID(),l);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (loadAllLabs != null) { try { loadAllLabs.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allLabsLoaded = true;
        return values;
    }	
	
    /**
     * Insert Lab by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertLab(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into lab (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into lab table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }

    
    //////////////////
    // CellLine stuff
    //////////////////

    /**
     * Load single CellLine by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public CellLine loadCellLine(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(cellNames.containsKey(name) && !forceDatabaseRefresh) { return cellNames.get(name); }
    	
    	CellLine c=null;
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from cellline where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
                c = new CellLine(rs);
                
                if(!cellIDs.containsKey(c.getDBID())) { 
                    cellIDs.put(c.getDBID(), c);
                    cellNames.put(c.getName(), c);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(c==null && insertIfNone){
			int id = insertCellLine(name);
			return loadCellLine(id, forceDatabaseRefresh);
		}else{
			return c;
		}
    }
    
    /**
     * Load single CellLine by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public CellLine loadCellLine(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(cellIDs.containsKey(dbid) && !forceDatabaseRefresh) { return cellIDs.get(dbid); }

        CellLine c = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from cellline where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            if(rs.next()) { 
                c = new CellLine(rs);
                cellIDs.put(dbid, c);
                cellNames.put(c.getName(), c);
            } else {
                throw new IllegalArgumentException("Unknown Cells DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return c;
    }
    
    /**
     * Load a collection of CellLines by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<CellLine> loadCellLines(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<CellLine> values = new LinkedList<CellLine>();
        for(int dbid : dbids) { values.addLast(loadCellLine(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all CellLines
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<CellLine> loadAllCellLines(boolean forceDatabaseRefresh) throws SQLException {
    	if(allCellsLoaded && !forceDatabaseRefresh){ return cellNames.values();}
    	if(forceDatabaseRefresh){cellNames.clear(); cellIDs.clear();}
        
    	HashSet<CellLine> values = new HashSet<CellLine>();
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from cellline");
	        rs = ps.executeQuery();
	
	        while(rs.next()) { 
	            CellLine c = new CellLine(rs);
	            values.add(c);
	            cellNames.put(c.getName(), c);
	            cellIDs.put(c.getDBID(),c);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allCellsLoaded = true;
        return values;
    }	

    /**
     * Insert a CellLine by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertCellLine(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into cellline (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into cellline table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }

    
    //////////////////
    // ExptCondition stuff
    //////////////////
    
    /**
     * Load a single ExptCondition by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ExptCondition loadExptCondition(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(condNames.containsKey(name) && !forceDatabaseRefresh) { return condNames.get(name); }
    	ExptCondition c=null;
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from exptcondition where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
                c = new ExptCondition(rs);
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(c==null && insertIfNone){
			int id = insertExptCondition(name);
			return loadExptCondition(id, forceDatabaseRefresh);
		}else{
			return c;
		}
    }        

    /**
     * Load a single ExptCondition by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ExptCondition loadExptCondition(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(condIDs.containsKey(dbid) && !forceDatabaseRefresh) {  return condIDs.get(dbid); }
        
        ExptCondition c = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from exptcondition where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
                c = new ExptCondition(rs);
                condIDs.put(dbid, c);
                condNames.put(c.getName(), c);
            } else {
                throw new IllegalArgumentException("Unknown ExptCondition DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return c;
    }

    /**
     * Load a collection of ExptConditions by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ExptCondition> loadExptConditions(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<ExptCondition> values = new LinkedList<ExptCondition>();
        for(int dbid : dbids) { values.addLast(loadExptCondition(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all ExptConditions
     * @param forceDatabaseRefresh
     * @return
     * @throws SQLException
     */
    public Collection<ExptCondition> loadAllExptConditions(boolean forceDatabaseRefresh) throws SQLException { 
    	if(allCondsLoaded && !forceDatabaseRefresh){ return condNames.values();}
    	if(forceDatabaseRefresh){condNames.clear(); condIDs.clear();}
        
    	HashSet<ExptCondition> values = new HashSet<ExptCondition>();
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from exptcondition");
	        rs = ps.executeQuery();
	
	        while(rs.next()) { 
	            ExptCondition c = new ExptCondition(rs);
	            values.add(c);
	            condNames.put(c.getName(), c);
	            condIDs.put(c.getDBID(),c);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allCondsLoaded = true;
        return values;
    }	

    /**
     * Insert ExptCondition by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertExptCondition(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into exptcondition (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into exptcondition table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }

    
    //////////////////
    // ExptTarget stuff
    //////////////////
    
    /**
     * Load a single ExptTarget by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ExptTarget loadExptTarget(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(targetNames.containsKey(name) && !forceDatabaseRefresh) { return targetNames.get(name); }
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		ExptTarget c=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from expttarget where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
                c = new ExptTarget(rs);
                if(!targetIDs.containsKey(c.getDBID())) { 
                    targetIDs.put(c.getDBID(), c);
                    targetNames.put(c.getName(), c);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(c==null && insertIfNone){
			int id = insertExptTarget(name);
			return loadExptTarget(id, forceDatabaseRefresh);
		}else{
			return c;
		}
    }
    
    /**
     * Load a single ExptTarget by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ExptTarget loadExptTarget(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(targetIDs.containsKey(dbid) && !forceDatabaseRefresh) { return targetIDs.get(dbid); }
        ExptTarget c = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from expttarget where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
                c = new ExptTarget(rs);
                targetIDs.put(dbid, c);
                targetNames.put(c.getName(), c);
            } else {
                throw new IllegalArgumentException("Unknown ExptTarget DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return c;
    }
  	
    /**
     * Load a collection of ExptTargets by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ExptTarget> loadExptTargets(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<ExptTarget> values = new LinkedList<ExptTarget>();
        for(int dbid : dbids) { values.addLast(loadExptTarget(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all ExptTargets
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ExptTarget> loadAllExptTargets(boolean forceDatabaseRefresh) throws SQLException { 
    	if(allTargetsLoaded && !forceDatabaseRefresh){ return targetNames.values();}
    	if(forceDatabaseRefresh){targetNames.clear(); targetIDs.clear();}
        
    	HashSet<ExptTarget> values = new HashSet<ExptTarget>();
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from expttarget");
	        rs = ps.executeQuery();
	
	        while(rs.next()) { 
	            ExptTarget f = new ExptTarget(rs);
	            values.add(f);
	            targetNames.put(f.getName(), f);
	            targetIDs.put(f.getDBID(),f);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allTargetsLoaded = true;
        return values;
    }	
	
    /**
     * Insert an ExptTarget by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertExptTarget(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into expttarget (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into expttarget table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }

    //////////////////
    // ExptType stuff
    //////////////////
    
    /**
     * Load a single ExptType by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ExptType loadExptType(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(exptTypeNames.containsKey(name) && !forceDatabaseRefresh) { return exptTypeNames.get(name); }
    	ExptType e=null;
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from expttype where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
            	e = new ExptType(rs);
                if(!exptTypeIDs.containsKey(e.getDBID())) { 
                    exptTypeIDs.put(e.getDBID(), e);
                    exptTypeNames.put(e.getName(), e);
                }
            }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(e==null && insertIfNone){
			int id = insertExptType(name);
			return loadExptType(id, forceDatabaseRefresh);
		}else{
			return e;
		}
    }
    
    /**
     * Load a single ExptType by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ExptType loadExptType(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(exptTypeIDs.containsKey(dbid) && !forceDatabaseRefresh) { return exptTypeIDs.get(dbid); }
        ExptType e = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from expttype where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            if(rs.next()) { 
                e = new ExptType(rs);
                exptTypeIDs.put(dbid, e);
                exptTypeNames.put(e.getName(), e);
            } else {
                throw new IllegalArgumentException("Unknown ExptType DBID: " + dbid);
            }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return e;
    }
    
    /**
     * Load a collection of ExptTypes by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ExptType> loadExptTypes(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<ExptType> values = new LinkedList<ExptType>();
        for(int dbid : dbids) { values.addLast(loadExptType(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all ExptTypes
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ExptType> loadAllExptTypes(boolean forceDatabaseRefresh) throws SQLException {
    	if(allExptTypesLoaded && !forceDatabaseRefresh){ return exptTypeNames.values();}
    	if(forceDatabaseRefresh){exptTypeNames.clear(); exptTypeIDs.clear();}
        
    	HashSet<ExptType> values = new HashSet<ExptType>();
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from expttype");
	        rs = ps.executeQuery();
	
	        while(rs.next()) { 
	        	ExptType e = new ExptType(rs);
	            values.add(e);
	            exptTypeNames.put(e.getName(), e);
	            exptTypeIDs.put(e.getDBID(),e);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allExptTypesLoaded = true;
        return values;
    }	
	
    /**
     * Insert an ExptType by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertExptType(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into expttype (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into expttype table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }

    //////////////////
    // ReadType stuff
    //////////////////
    
    /**
     * Load a single ReadType by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ReadType loadReadType(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException {
    	if(readTypeNames.containsKey(name) && !forceDatabaseRefresh) { return readTypeNames.get(name); }
    	ReadType r =null;
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from readtype where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
            	r = new ReadType(rs);
                if(!readTypeIDs.containsKey(r.getDBID())) { 
                    readTypeIDs.put(r.getDBID(), r);
                    readTypeNames.put(r.getName(), r);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(r==null && insertIfNone){
			int id = insertReadType(name);
			return loadReadType(id, forceDatabaseRefresh);
		}else{
			return r;
		}
    }
        
    /**
     * Load a single ReadType by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public ReadType loadReadType(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(readTypeIDs.containsKey(dbid) && !forceDatabaseRefresh) { return readTypeIDs.get(dbid); }

        ReadType e = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from readtype where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            if(rs.next()) { 
                e = new ReadType(rs);
                readTypeIDs.put(dbid, e);
                readTypeNames.put(e.getName(), e);
            } else {
                throw new IllegalArgumentException("Unknown ReadType DBID: " + dbid);
            }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }       
        return e;
    }
    
    /**
     * Load a collection of ReadTypes by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ReadType> loadReadTypes(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<ReadType> values = new LinkedList<ReadType>();
        for(int dbid : dbids) { values.addLast(loadReadType(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all ReadTypes
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<ReadType> loadAllReadTypes(boolean forceDatabaseRefresh) throws SQLException {
    	if(allReadTypesLoaded && !forceDatabaseRefresh){ return readTypeNames.values();}
    	if(forceDatabaseRefresh){readTypeNames.clear(); readTypeIDs.clear();}
    	
    	HashSet<ReadType> values = new HashSet<ReadType>();
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from readtype");
	        rs = ps.executeQuery();
	
	        while(rs.next()) { 
	        	ReadType r = new ReadType(rs);
	            values.add(r);
	            readTypeNames.put(r.getName(), r);
	            readTypeIDs.put(r.getDBID(),r);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allReadTypesLoaded = true;
        return values;
    }	
	
    /**
     * Insert a ReadType by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertReadType(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into readtype (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into readtype table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }


    //////////////////
    // AlignType stuff
    //////////////////
    
    /**
     * Load a single AlignType by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public AlignType loadAlignType(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(alignTypeNames.containsKey(name) && !forceDatabaseRefresh) { return alignTypeNames.get(name); }
    	
    	AlignType a =null;
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from aligntype where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
            	a = new AlignType(rs);
                
                if(!alignTypeIDs.containsKey(a.getDBID())) { 
                    alignTypeIDs.put(a.getDBID(), a);
                    alignTypeNames.put(a.getName(), a);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(a==null && insertIfNone){
			int id = insertAlignType(name);
			return loadAlignType(id, forceDatabaseRefresh);
		}else{
			return a;
		}
    }
    
    /**
     * Load a single AlignType by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public AlignType loadAlignType(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
        if(alignTypeIDs.containsKey(dbid) && !forceDatabaseRefresh) { return alignTypeIDs.get(dbid); }
        
        AlignType a = null;
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from aligntype where id=?");
            ps.setInt(1, dbid);
            rs = ps.executeQuery();
            if(rs.next()) { 
                a = new AlignType(rs);
                alignTypeIDs.put(dbid, a);
                alignTypeNames.put(a.getName(), a);
            } else {
                throw new IllegalArgumentException("Unknown AlignType DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return a;
    }
    
    /**
     * Load a collection of AlignTypes by IDs
     * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<AlignType> loadAlignTypes(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
        LinkedList<AlignType> values = new LinkedList<AlignType>();
        for(int dbid : dbids) { values.addLast(loadAlignType(dbid, forceDatabaseRefresh)); }
        return values;
    }

    /**
     * Load all AlignTypes
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public Collection<AlignType> loadAllAlignTypes(boolean forceDatabaseRefresh) throws SQLException {
    	if(allAlignTypesLoaded && !forceDatabaseRefresh){ return alignTypeNames.values();}
    	if(forceDatabaseRefresh){alignTypeNames.clear(); alignTypeIDs.clear();}
    	
    	HashSet<AlignType> values = new HashSet<AlignType>();
        Connection cxn = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name from aligntype");
	        rs = ps.executeQuery();
	
	        while(rs.next()) { 
	        	AlignType a = new AlignType(rs);
	            values.add(a);
	            alignTypeNames.put(a.getName(), a);
	            alignTypeIDs.put(a.getDBID(),a);
	        }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allAlignTypesLoaded = true;
        return values;
    }	
	
    /**
     * Insert an AlignType by name
     * @param n
     * @return
     * @throws SQLException
     */
    private int insertAlignType(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
	    						    java.sql.ResultSet.CONCUR_UPDATABLE);
	        s.executeUpdate("insert into aligntype (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
	        rs = s.getGeneratedKeys();

	        if (rs.next())
	            id = rs.getInt(1);
	        else 
	        	throw new IllegalArgumentException("Unable to insert new entry into aligntype table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
        return id;
    }

	//////////////////
	// SeqDataUser stuff
	//////////////////
    
    /**
     * Load a single SeqDataUser by name
     * @param name
     * @param insertIfNone : insert into the database if there is no such entry
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
    public SeqDataUser loadSeqDataUser(String name, boolean insertIfNone, boolean forceDatabaseRefresh) throws SQLException { 
    	if(seqDataUserNames.containsKey(name) && !forceDatabaseRefresh) { return seqDataUserNames.get(name); }
    	SeqDataUser a=null;
    	Connection cxn = null;
    	PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name, admin from seqdatauser where name=?");
            ps.setString(1, name);
            rs = ps.executeQuery();
            
            if(rs.next()) { 
            	a = new SeqDataUser(rs);
                if(!seqDataUserIDs.containsKey(a.getDBID())) { 
                    seqDataUserIDs.put(a.getDBID(), a);
                    seqDataUserNames.put(a.getName(), a);
                }
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(a==null && insertIfNone){
			int id = insertSeqDataUser(name);
			return loadSeqDataUser(id, forceDatabaseRefresh);
		}else{
			return a;
		}
    }
    
    /**
     * Load a single SeqDataUser by ID
     * @param dbid
     * @param forceDatabaseRefresh : ignore cache and query database directly
     * @return
     * @throws SQLException
     */
	public SeqDataUser loadSeqDataUser(int dbid, boolean forceDatabaseRefresh) throws SQLException { 
		if(seqDataUserIDs.containsKey(dbid) && !forceDatabaseRefresh) { return seqDataUserIDs.get(dbid); }
	
		SeqDataUser a = null;
		Connection cxn = null;
		PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name, admin from seqdatauser where id=?");
			ps.setInt(1, dbid);
			rs = ps.executeQuery();
			if(rs.next()) { 
				a = new SeqDataUser(rs);
				seqDataUserIDs.put(dbid, a);
				seqDataUserNames.put(a.getName(), a);
			} else {
				throw new IllegalArgumentException("Unknown AlignType DBID: " + dbid);
			}
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }    
		return a;
	}
	
	/**
	 * Load a collection of SeqDataUsers by IDs
	 * @param dbids
     * @param forceDatabaseRefresh : ignore cache and query database directly
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqDataUser> loadSeqDataUsers(Collection<Integer> dbids, boolean forceDatabaseRefresh) throws SQLException {
		LinkedList<SeqDataUser> values = new LinkedList<SeqDataUser>();
		for(int dbid : dbids) { values.addLast(loadSeqDataUser(dbid, forceDatabaseRefresh)); }
		return values;
	}
	
	/**
	 * Load all SeqDataUsers
     * @param forceDatabaseRefresh : ignore cache and query database directly
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqDataUser> loadAllSeqDataUsers(boolean forceDatabaseRefresh) throws SQLException {
		if(allSeqDataUsersLoaded && !forceDatabaseRefresh){ return seqDataUserNames.values();}
		if(forceDatabaseRefresh){seqDataUserNames.clear(); seqDataUserIDs.clear();}
		
		HashSet<SeqDataUser> values = new HashSet<SeqDataUser>();
		Connection cxn = null;
		PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = cxn.prepareStatement("select id, name, admin from seqdatauser");
			rs = ps.executeQuery();
		
			while(rs.next()) { 
				SeqDataUser a = new SeqDataUser(rs);
				values.add(a);
				seqDataUserNames.put(a.getName(), a);
				seqDataUserIDs.put(a.getDBID(),a);
			}
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		allSeqDataUsersLoaded = true;
		return values;
	}	
	
	/**
	 * Insert a SeqDataUser by name
	 * @param n
	 * @return
	 * @throws SQLException
	 */
	private int insertSeqDataUser(String n) throws SQLException {
		Statement s = null;
		ResultSet rs = null;
		int id=-1;
		Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			s = cxn.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY,
						    	java.sql.ResultSet.CONCUR_UPDATABLE);
			s.executeUpdate("insert into seqdatauser (name) values ('" + n + "')", Statement.RETURN_GENERATED_KEYS);
			rs = s.getGeneratedKeys();
	
			if (rs.next())
				id = rs.getInt(1);
			else 
				throw new IllegalArgumentException("Unable to insert new entry into seqdatauser table"); 
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (s != null) { try { s.close();} catch (SQLException ex) { } }
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
	    }
		return id;
	}

}
