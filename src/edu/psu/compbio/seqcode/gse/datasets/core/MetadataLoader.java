package edu.psu.compbio.seqcode.gse.datasets.core;

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

import edu.psu.compbio.seqcode.gse.utils.database.DatabaseConnectionManager;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/**
 * MetadataLoader is the interface for interacting with metadata entries in the core database.
 * 
 * 
 * To make MetadataLoader thread-safe, many of its method have an internal block
   that is synchronized on some PreparedStatement.  To prevent deadlocks, we want a 
   directed acyclic graph of which Loader uses which.

   Please add entries below for other loaders that also synchronize on themselves or fields
   to ensure that we don't create cycles (indentation below an entry indicates "uses").

   MetadataLoader
   ExpressionMetadataLoader
   ExpressionLoader
       ProbeMappingLoader
   TimeSeriesLoader
       ExpressionLoader
           ProbeMappingLoader
   WeightMatrixLoader
   SeqDataLoader
   
*/

public class MetadataLoader{
    
    public static final String role = "core";

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
    	
    public MetadataLoader() throws SQLException { 
        
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
    }
	
    
 
    //////////////////
    // Lab stuff
    //////////////////
    
    public Lab getLab(String name) throws SQLException { 
    	Connection cxn = null;
    	Lab l = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadLabsByName = cxn.prepareStatement("select id, name from lab where name=?");
            loadLabsByName.setString(1, name);
            ResultSet rs = loadLabsByName.executeQuery();
            
            if(rs.next()) { 
                l = new Lab(rs);
                
                if(!labIDs.containsKey(l.getDBID())) { 
                    labIDs.put(l.getDBID(), l);
                    labNames.put(l.getName(), l);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(l == null){
			int id = insertLab(name);
			return loadLab(id);
		}else{
			return l;
		}
    }
    
    public Lab findLab(String name) throws SQLException { 
    	Connection cxn = null;
    	Lab l = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadLabsByName = cxn.prepareStatement("select id, name from lab where name=?");
            loadLabsByName.setString(1, name);
            ResultSet rs = loadLabsByName.executeQuery();
            
            if(rs.next()) { 
                l = new Lab(rs);
                
                if(!labIDs.containsKey(l.getDBID())) { 
                    labIDs.put(l.getDBID(), l);
                    labNames.put(l.getName(), l);
                }
            }            
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return l;
    }
    
    public Lab loadLab(int dbid) throws SQLException { 
        if(labIDs.containsKey(dbid)) { return labIDs.get(dbid); }

        Lab l = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadLabs = cxn.prepareStatement("select id, name from lab where id=?");
            loadLabs.setInt(1, dbid);
            ResultSet rs = loadLabs.executeQuery();
            if(rs.next()) { 
                l = new Lab(rs);
                labIDs.put(dbid, l);
                labNames.put(l.getName(), l);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Lab DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return l;
    }
    
    public Collection<Lab> loadAllLabs(Collection<Integer> dbids) throws SQLException {

        LinkedList<Lab> values = new LinkedList<Lab>();
        for(int dbid : dbids) { values.addLast(loadLab(dbid)); }
        return values;
    }

    public Collection<Lab> loadAllLabs() throws SQLException {
        HashSet<Lab> values = new HashSet<Lab>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllLabs = cxn.prepareStatement("select id, name from lab");
	        ResultSet rs = loadAllLabs.executeQuery();
	
	        while(rs.next()) { 
	        	Lab l = new Lab(rs);
	            values.add(l);
	            labNames.put(l.getName(), l);
	            labIDs.put(l.getDBID(),l);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	
	
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
	        rs.close();
	        rs = null;
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

    public CellLine getCellLine(String name) throws SQLException { 
    	CellLine c=null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadCellsByName = cxn.prepareStatement("select id, name from cellline where name=?");
            loadCellsByName.setString(1, name);
            ResultSet rs = loadCellsByName.executeQuery();
            
            if(rs.next()) { 
                c = new CellLine(rs);
                
                if(!cellIDs.containsKey(c.getDBID())) { 
                    cellIDs.put(c.getDBID(), c);
                    cellNames.put(c.getName(), c);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(c==null){
			int id = insertCellLine(name);
			return loadCellLine(id);
		}else{
			return c;
		}
    }
    
    public CellLine findCellLine(String name) throws SQLException { 
    	CellLine c = null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadCellsByName = cxn.prepareStatement("select id, name from cellline where name=?");
            loadCellsByName.setString(1, name);
            ResultSet rs = loadCellsByName.executeQuery();
            
            if(rs.next()) { 
                c = new CellLine(rs);
                if(!cellIDs.containsKey(c.getDBID())) { 
                    cellIDs.put(c.getDBID(), c);
                    cellNames.put(c.getName(), c);
                }
            }            
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return c;
    }
    
    public CellLine loadCellLine(int dbid) throws SQLException { 
        if(cellIDs.containsKey(dbid)) { return cellIDs.get(dbid); }

        CellLine c = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadCells = cxn.prepareStatement("select id, name from cellline where id=?");
            loadCells.setInt(1, dbid);
            ResultSet rs = loadCells.executeQuery();
            if(rs.next()) { 
                c = new CellLine(rs);
                cellIDs.put(dbid, c);
                cellNames.put(c.getName(), c);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Cells DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return c;
    }
    
    public Collection<CellLine> loadAllCellLines(Collection<Integer> dbids) throws SQLException {

        LinkedList<CellLine> values = new LinkedList<CellLine>();
        for(int dbid : dbids) { values.addLast(loadCellLine(dbid)); }
        return values;
    }

    public Collection<CellLine> loadAllCellLines() throws SQLException {
        HashSet<CellLine> values = new HashSet<CellLine>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllCells = cxn.prepareStatement("select id, name from cellline");
	        ResultSet rs = loadAllCells.executeQuery();
	
	        while(rs.next()) { 
	            CellLine c = new CellLine(rs);
	            values.add(c);
	            cellNames.put(c.getName(), c);
	            cellIDs.put(c.getDBID(),c);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	

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
	        rs.close();
	        rs = null;
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
    
    public ExptCondition getExptCondition(String name) throws SQLException { 
    	ExptCondition c=null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadCondByName = cxn.prepareStatement("select id, name from exptcondition where name=?");
            loadCondByName.setString(1, name);
            ResultSet rs = loadCondByName.executeQuery();
            
            if(rs.next()) { 
                c = new ExptCondition(rs);
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(c==null){
			int id = insertExptCondition(name);
			return loadExptCondition(id);
		}else{
			return c;
		}
    }        

    public ExptCondition findExptCondition(String name) throws SQLException { 
    	ExptCondition c=null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadCondByName = cxn.prepareStatement("select id, name from exptcondition where name=?");
            loadCondByName.setString(1, name);
            ResultSet rs = loadCondByName.executeQuery();
            
            if(rs.next()) { 
                c = new ExptCondition(rs);
                
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
            }            
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return c;
    }
    
    public ExptCondition loadExptCondition(int dbid) throws SQLException { 
        if(condIDs.containsKey(dbid)) {  return condIDs.get(dbid); }
        
        ExptCondition c = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadCond = cxn.prepareStatement("select id, name from exptcondition where id=?");
            loadCond.setInt(1, dbid);
            ResultSet rs = loadCond.executeQuery();
            
            if(rs.next()) { 
                c = new ExptCondition(rs);
                condIDs.put(dbid, c);
                condNames.put(c.getName(), c);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ExptCondition DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return c;
    }

    public Collection<ExptCondition> loadAllExptConditions(Collection<Integer> dbids) throws SQLException {

        LinkedList<ExptCondition> values = new LinkedList<ExptCondition>();
        for(int dbid : dbids) { values.addLast(loadExptCondition(dbid)); }
        return values;
    }

    public Collection<ExptCondition> loadAllExptConditions() throws SQLException { 
        HashSet<ExptCondition> values = new HashSet<ExptCondition>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllCond = cxn.prepareStatement("select id, name from exptcondition");
	        ResultSet rs = loadAllCond.executeQuery();
	
	        while(rs.next()) { 
	            ExptCondition c = new ExptCondition(rs);
	            values.add(c);
	            condNames.put(c.getName(), c);
	            condIDs.put(c.getDBID(),c);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	

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
	        rs.close();
	        rs = null;
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
    
    public ExptTarget getExptTarget(String name) throws SQLException { 
    	Connection cxn = null;
    	ExptTarget c=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadTargetsByName = cxn.prepareStatement("select id, name from expttarget where name=?");
            loadTargetsByName.setString(1, name);
            ResultSet rs = loadTargetsByName.executeQuery();
            
            if(rs.next()) { 
                c = new ExptTarget(rs);
                if(!targetIDs.containsKey(c.getDBID())) { 
                    targetIDs.put(c.getDBID(), c);
                    targetNames.put(c.getName(), c);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(c==null){
			int id = insertFactor(name);
			return loadExptTarget(id);
		}else{
			return c;
		}
    }
    
    public ExptTarget findExptTarget(String name) throws SQLException { 
    	ExptTarget c = null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadTargetsByName = cxn.prepareStatement("select id, name from expttarget where name=?");
            loadTargetsByName.setString(1, name);
            ResultSet rs = loadTargetsByName.executeQuery();
            
            if(rs.next()) { 
                c = new ExptTarget(rs);
                if(!targetIDs.containsKey(c.getDBID())) { 
                    targetIDs.put(c.getDBID(), c);
                    targetNames.put(c.getName(), c);
                }
            }        
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return c;
    }
    
    public ExptTarget loadExptTarget(int dbid) throws SQLException { 
        if(targetIDs.containsKey(dbid)) { return targetIDs.get(dbid); }
        ExptTarget c = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadTargets = cxn.prepareStatement("select id, name from expttarget where id=?");
            loadTargets.setInt(1, dbid);
            ResultSet rs = loadTargets.executeQuery();
            
            if(rs.next()) { 
                c = new ExptTarget(rs);
                targetIDs.put(dbid, c);
                targetNames.put(c.getName(), c);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ExptTarget DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return c;
    }
  	
    public Collection<ExptTarget> loadAllExptTargets(Collection<Integer> dbids) throws SQLException {

        LinkedList<ExptTarget> values = new LinkedList<ExptTarget>();
        for(int dbid : dbids) { values.addLast(loadExptTarget(dbid)); }
        return values;
    }

    public Collection<ExptTarget> loadAllExptTargets() throws SQLException { 
        HashSet<ExptTarget> values = new HashSet<ExptTarget>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllTargets = cxn.prepareStatement("select id, name from expttarget");
	        ResultSet rs = loadAllTargets.executeQuery();
	
	        while(rs.next()) { 
	            ExptTarget f = new ExptTarget(rs);
	            values.add(f);
	            targetNames.put(f.getName(), f);
	            targetIDs.put(f.getDBID(),f);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	
	
    private int insertFactor(String n) throws SQLException {
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
	        rs.close();
	        rs = null;
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
    
    public ExptType getExptType(String name) throws SQLException { 
    	ExptType e=null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadExptTypesByName = cxn.prepareStatement("select id, name from expttype where name=?");
            loadExptTypesByName.setString(1, name);
            ResultSet rs = loadExptTypesByName.executeQuery();
            
            if(rs.next()) { 
            	e = new ExptType(rs);
                if(!exptTypeIDs.containsKey(e.getDBID())) { 
                    exptTypeIDs.put(e.getDBID(), e);
                    exptTypeNames.put(e.getName(), e);
                }
            }
            rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(e==null){
			int id = insertExptType(name);
			return loadExptType(id);
		}else{
			return e;
		}
    }
    
    public ExptType findExptType(String name) throws SQLException { 
    	ExptType e = null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadExptTypesByName = cxn.prepareStatement("select id, name from expttype where name=?");
            loadExptTypesByName.setString(1, name);
            ResultSet rs = loadExptTypesByName.executeQuery();
            
            if(rs.next()) { 
            	e = new ExptType(rs);
                if(!exptTypeIDs.containsKey(e.getDBID())) { 
                    exptTypeIDs.put(e.getDBID(), e);
                    exptTypeNames.put(e.getName(), e);
                }
            }            
            rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return e;
    }
    
    public ExptType loadExptType(int dbid) throws SQLException { 
        if(exptTypeIDs.containsKey(dbid)) { return exptTypeIDs.get(dbid); }
        ExptType e = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadExptTypes = cxn.prepareStatement("select id, name from expttype where id=?");
            loadExptTypes.setInt(1, dbid);
            ResultSet rs = loadExptTypes.executeQuery();
            if(rs.next()) { 
                e = new ExptType(rs);
                exptTypeIDs.put(dbid, e);
                exptTypeNames.put(e.getName(), e);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ExptType DBID: " + dbid);
            }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return e;
    }
    
    public Collection<ExptType> loadAllExptTypes(Collection<Integer> dbids) throws SQLException {

        LinkedList<ExptType> values = new LinkedList<ExptType>();
        for(int dbid : dbids) { values.addLast(loadExptType(dbid)); }
        return values;
    }

    public Collection<ExptType> loadAllExptTypes() throws SQLException {
        HashSet<ExptType> values = new HashSet<ExptType>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllExptTypes = cxn.prepareStatement("select id, name from expttype");
	        ResultSet rs = loadAllExptTypes.executeQuery();
	
	        while(rs.next()) { 
	        	ExptType e = new ExptType(rs);
	            values.add(e);
	            exptTypeNames.put(e.getName(), e);
	            exptTypeIDs.put(e.getDBID(),e);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	
	
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
	        rs.close();
	        rs = null;
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
    
    public ReadType getReadType(String name) throws SQLException { 
    	ReadType r =null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadReadTypesByName = cxn.prepareStatement("select id, name from readtype where name=?");
            loadReadTypesByName.setString(1, name);
            ResultSet rs = loadReadTypesByName.executeQuery();
            
            if(rs.next()) { 
            	r = new ReadType(rs);
                if(!readTypeIDs.containsKey(r.getDBID())) { 
                    readTypeIDs.put(r.getDBID(), r);
                    readTypeNames.put(r.getName(), r);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(r==null){
			int id = insertReadType(name);
			return loadReadType(id);
		}else{
			return r;
		}
    }
    
    public ReadType findReadType(String name) throws SQLException { 
    	ReadType r = null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadReadTypesByName = cxn.prepareStatement("select id, name from readtype where name=?");
            loadReadTypesByName.setString(1, name);
            ResultSet rs = loadReadTypesByName.executeQuery();
            
            if(rs.next()) { 
            	r = new ReadType(rs);
                if(!readTypeIDs.containsKey(r.getDBID())) { 
                    readTypeIDs.put(r.getDBID(), r);
                    readTypeNames.put(r.getName(), r);
                }
            }            
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return r;
    }
    
    public ReadType loadReadType(int dbid) throws SQLException { 
        if(readTypeIDs.containsKey(dbid)) { return readTypeIDs.get(dbid); }

        ReadType e = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadReadTypes = cxn.prepareStatement("select id, name from readtype where id=?");
            loadReadTypes.setInt(1, dbid);
            ResultSet rs = loadReadTypes.executeQuery();
            if(rs.next()) { 
                e = new ReadType(rs);
                readTypeIDs.put(dbid, e);
                readTypeNames.put(e.getName(), e);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ReadType DBID: " + dbid);
            }
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }       
        return e;
    }
    
    public Collection<ReadType> loadAllReadTypes(Collection<Integer> dbids) throws SQLException {

        LinkedList<ReadType> values = new LinkedList<ReadType>();
        for(int dbid : dbids) { values.addLast(loadReadType(dbid)); }
        return values;
    }

    public Collection<ReadType> loadAllReadTypes() throws SQLException {
        HashSet<ReadType> values = new HashSet<ReadType>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllReadTypes = cxn.prepareStatement("select id, name from readtype");
	        ResultSet rs = loadAllReadTypes.executeQuery();
	
	        while(rs.next()) { 
	        	ReadType r = new ReadType(rs);
	            values.add(r);
	            readTypeNames.put(r.getName(), r);
	            readTypeIDs.put(r.getDBID(),r);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	
	
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
	        rs.close();
	        rs = null;
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
    
    public AlignType getAlignType(String name) throws SQLException { 
    	AlignType a =null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAlignTypesByName = cxn.prepareStatement("select id, name from aligntype where name=?");
            loadAlignTypesByName.setString(1, name);
            ResultSet rs = loadAlignTypesByName.executeQuery();
            
            if(rs.next()) { 
            	a = new AlignType(rs);
                
                if(!alignTypeIDs.containsKey(a.getDBID())) { 
                    alignTypeIDs.put(a.getDBID(), a);
                    alignTypeNames.put(a.getName(), a);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(a==null){
			int id = insertAlignType(name);
			return loadAlignType(id);
		}else{
			return a;
		}
    }
    
    public AlignType findAlignType(String name) throws SQLException { 
    	AlignType a = null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAlignTypesByName = cxn.prepareStatement("select id, name from aligntype where name=?");
            loadAlignTypesByName.setString(1, name);
            ResultSet rs = loadAlignTypesByName.executeQuery();
            
            if(rs.next()) { 
            	a = new AlignType(rs);
                
                if(!alignTypeIDs.containsKey(a.getDBID())) { 
                    alignTypeIDs.put(a.getDBID(), a);
                    alignTypeNames.put(a.getName(), a);
                }
            }            
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return a;
    }
    
    public AlignType loadAlignType(int dbid) throws SQLException { 
        if(alignTypeIDs.containsKey(dbid)) { return alignTypeIDs.get(dbid); }
        AlignType a = null;
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAlignTypes = cxn.prepareStatement("select id, name from aligntype where id=?");
            loadAlignTypes.setInt(1, dbid);
            ResultSet rs = loadAlignTypes.executeQuery();
            if(rs.next()) { 
                a = new AlignType(rs);
                alignTypeIDs.put(dbid, a);
                alignTypeNames.put(a.getName(), a);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown AlignType DBID: " + dbid);
            }
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return a;
    }
    
    public Collection<AlignType> loadAllAlignTypes(Collection<Integer> dbids) throws SQLException {

        LinkedList<AlignType> values = new LinkedList<AlignType>();
        for(int dbid : dbids) { values.addLast(loadAlignType(dbid)); }
        return values;
    }

    public Collection<AlignType> loadAllAlignTypes() throws SQLException {
        HashSet<AlignType> values = new HashSet<AlignType>();
        Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllAlignTypes = cxn.prepareStatement("select id, name from aligntype");
	        ResultSet rs = loadAllAlignTypes.executeQuery();
	
	        while(rs.next()) { 
	        	AlignType a = new AlignType(rs);
	            values.add(a);
	            alignTypeNames.put(a.getName(), a);
	            alignTypeIDs.put(a.getDBID(),a);
	        }
	        rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
        return values;
    }	
	
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
	        rs.close();
	        rs = null;
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
    public SeqDataUser getSeqDataUser(String name) throws SQLException { 
    	SeqDataUser a=null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadSeqDataUsersByName = cxn.prepareStatement("select id, name, admin from seqdatauser where name=?");
            loadSeqDataUsersByName.setString(1, name);
            ResultSet rs = loadSeqDataUsersByName.executeQuery();
            
            if(rs.next()) { 
            	a = new SeqDataUser(rs);
                if(!seqDataUserIDs.containsKey(a.getDBID())) { 
                    seqDataUserIDs.put(a.getDBID(), a);
                    seqDataUserNames.put(a.getName(), a);
                }
            }
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		if(a==null){
			int id = insertSeqDataUser(name);
			return loadSeqDataUser(id);
		}else{
			return a;
		}
    }
    
    public SeqDataUser findSeqDataUser(String name) throws SQLException { 
    	SeqDataUser a = null;
    	Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadSeqDataUsersByName = cxn.prepareStatement("select id, name, admin from seqdatauser where name=?");
            loadSeqDataUsersByName.setString(1, name);
            ResultSet rs = loadSeqDataUsersByName.executeQuery();
            
            if(rs.next()) { 
            	a = new SeqDataUser(rs);
                if(!seqDataUserIDs.containsKey(a.getDBID())) { 
                    seqDataUserIDs.put(a.getDBID(), a);
                    seqDataUserNames.put(a.getName(), a);
                }
            }            
            rs.close();
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return a;
    }
	public SeqDataUser loadSeqDataUser(int dbid) throws SQLException { 
		if(seqDataUserIDs.containsKey(dbid)) { return seqDataUserIDs.get(dbid); }
	
		SeqDataUser a = null;
		Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadSeqDataUsers = cxn.prepareStatement("select id, name, admin from seqdatauser where id=?");
			loadSeqDataUsers.setInt(1, dbid);
			ResultSet rs = loadSeqDataUsers.executeQuery();
			if(rs.next()) { 
				a = new SeqDataUser(rs);
				seqDataUserIDs.put(dbid, a);
				seqDataUserNames.put(a.getName(), a);
				rs.close();
			} else {
				rs.close();
				throw new IllegalArgumentException("Unknown AlignType DBID: " + dbid);
			}
		} catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }    
		return a;
	}
	
	public Collection<SeqDataUser> loadAllSeqDataUsers(Collection<Integer> dbids) throws SQLException {
		LinkedList<SeqDataUser> values = new LinkedList<SeqDataUser>();
		for(int dbid : dbids) { values.addLast(loadSeqDataUser(dbid)); }
		return values;
	}
	
	public Collection<SeqDataUser> loadAllSeqDataUsers() throws SQLException {
		HashSet<SeqDataUser> values = new HashSet<SeqDataUser>();
		Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement loadAllSeqDataUsers = cxn.prepareStatement("select id, name, admin from seqdatauser");
			ResultSet rs = loadAllSeqDataUsers.executeQuery();
		
			while(rs.next()) { 
				SeqDataUser a = new SeqDataUser(rs);
				values.add(a);
				seqDataUserNames.put(a.getName(), a);
				seqDataUserIDs.put(a.getDBID(),a);
			}
			rs.close();
		} catch (UnknownRoleException ex) {
            throw new IllegalArgumentException("Unknown role: " + role, ex);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		return values;
	}	
	
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
			rs.close();
			rs = null;
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
