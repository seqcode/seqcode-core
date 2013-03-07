package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/* to make MetadataLoader thread-safe, many of its method have an internal block
   that is synchronized on some PreparedStatement.  To prevent deadlocks, we want a 
   directed acyclic graph of which Loader uses which.

   Please add entries below for other loaders that also synchronize on themselves or fields
   to ensure that we don't create cycles (indentation below an entry indicates "uses").

   MetadataLoader
       ChipChipMetadataLoader
       ExpressionMetadataLoader
       ExpressionLoader
           ProbeMappingLoader
   TimeSeriesLoader
       ExpressionLoader
           ProbeMappingLoader
   ChipPetLoader
   BindingScanLoader
   LocatorLoader
   OrthologyLoader (uses itself through OrthologyPair constructor)
   TextFunctionLoader
   function.DatabaseFunctionLoader
   WeightMatrixLoader
   DistributionLoader
   AlignmentLoader
       ConservationLoader (uses an AlignmentLoader)
   

*/

public class MetadataLoader implements edu.psu.compbio.seqcode.gse.utils.Closeable {
    
    public static final String role = "core";

    private Map<String,Lab> labNames;
    private Map<String,CellLine> cellNames;
    private Map<String,ExptCondition> condNames;
    private Map<String,ExptTarget> targetNames;
    private Map<String,ExptType> exptTypeNames;
    private Map<String,ReadType> readTypeNames;
    private Map<String,AlignType> alignTypeNames;
	
    private Map<Integer,Lab> labIDs;
    private Map<Integer,CellLine> cellIDs;
    private Map<Integer,ExptCondition> condIDs;
    private Map<Integer,ExptTarget> targetIDs;
    private Map<Integer,ExptType> exptTypeIDs;
    private Map<Integer,ReadType> readTypeIDs;
    private Map<Integer,AlignType> alignTypeIDs;
	
    private java.sql.Connection cxn;
    
    private PreparedStatement loadLabs, loadCells, loadCond, loadTargets, loadExptTypes, loadReadTypes, loadAlignTypes;
    private PreparedStatement loadAllLabs, loadAllCells, loadAllCond, loadAllTargets, loadAllExptTypes, loadAllReadTypes, loadAllAlignTypes;
    private PreparedStatement loadLabsByName, loadCellsByName, loadCondByName, loadTargetsByName, loadExptTypesByName, loadReadTypesByName, loadAlignTypesByName;
	
    public MetadataLoader() throws SQLException { 
        try {
            cxn = DatabaseFactory.getConnection(role);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }
        
        loadLabs = cxn.prepareStatement("select id, name from lab where id=?");
        loadCells = cxn.prepareStatement("select id, name from cellline where id=?");
        loadCond = cxn.prepareStatement("select id, name from exptcondition where id=?");
        loadTargets = cxn.prepareStatement("select id, name from expttarget where id=?");
        loadExptTypes = cxn.prepareStatement("select id, name from expttype where id=?");
        loadReadTypes = cxn.prepareStatement("select id, name from readtype where id=?");
        loadAlignTypes = cxn.prepareStatement("select id, name from aligntype where id=?");

        loadAllLabs = cxn.prepareStatement("select id, name from lab");
        loadAllCells = cxn.prepareStatement("select id, name from cellline");
        loadAllCond = cxn.prepareStatement("select id, name from exptcondition");
        loadAllTargets = cxn.prepareStatement("select id, name from expttarget");
        loadAllExptTypes = cxn.prepareStatement("select id, name from expttype");
        loadAllReadTypes = cxn.prepareStatement("select id, name from readtype");
        loadAllAlignTypes = cxn.prepareStatement("select id, name from aligntype");

        loadLabsByName = cxn.prepareStatement("select id, name from lab where name=?");
        loadCellsByName = cxn.prepareStatement("select id, name from cellline where name=?");
        loadCondByName = cxn.prepareStatement("select id, name from exptcondition where name=?");
        loadTargetsByName = cxn.prepareStatement("select id, name from expttarget where name=?");
        loadExptTypesByName = cxn.prepareStatement("select id, name from expttype where name=?");
        loadReadTypesByName = cxn.prepareStatement("select id, name from readtype where name=?");
        loadAlignTypesByName = cxn.prepareStatement("select id, name from aligntype where name=?");

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
    }
	
    public boolean isClosed() { return cxn==null; }

    public void close() { 
        try {
        	loadLabs.close(); loadLabs=null;
            loadCells.close(); loadCells = null;
            loadCond.close();  loadCond = null;
            loadTargets.close(); loadTargets = null;
            loadExptTypes.close(); loadExptTypes = null;
            loadReadTypes.close(); loadReadTypes = null;
            loadAlignTypes.close(); loadAlignTypes = null;

            loadAllLabs.close(); loadAllLabs = null;
            loadAllCells.close(); loadAllCells = null;
            loadAllCond.close();  loadAllCond = null;
            loadAllTargets.close(); loadAllTargets = null;
            loadAllExptTypes.close(); loadAllExptTypes = null;
            loadAllReadTypes.close(); loadAllReadTypes = null;
            loadAllAlignTypes.close(); loadAllAlignTypes = null;
            
            loadLabsByName.close(); loadLabsByName=null;
            loadCellsByName.close();  loadCellsByName = null;
            loadCondByName.close(); loadCondByName = null;
            loadTargetsByName.close(); loadTargetsByName = null;
            loadExptTypesByName.close(); loadExptTypesByName = null;
            loadReadTypesByName.close(); loadReadTypesByName = null;
            loadAlignTypesByName.close(); loadAlignTypesByName = null;
        } catch (SQLException e) {
            e.printStackTrace();
        }
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }
    
    public java.sql.Connection getConnection() { return cxn; }
 
    //////////////////
    // Lab stuff
    //////////////////
    
    public Lab getLab(String name) throws SQLException { 
        synchronized (loadLabsByName) {
            loadLabsByName.setString(1, name);
            ResultSet rs = loadLabsByName.executeQuery();
            
            if(rs.next()) { 
                Lab l = new Lab(rs);
                rs.close();
                
                if(!labIDs.containsKey(l.getDBID())) { 
                    labIDs.put(l.getDBID(), l);
                    labNames.put(l.getName(), l);
                }
                
                return l;
            }
            rs.close();
        }
        int id = insertLab(name);
        return loadLab(id);
    }
    
    public Lab findLab(String name) throws SQLException { 
        synchronized (loadLabsByName) {
            loadLabsByName.setString(1, name);
            ResultSet rs = loadLabsByName.executeQuery();
            
            if(rs.next()) { 
                Lab l = new Lab(rs);
                rs.close();
                
                if(!labIDs.containsKey(l.getDBID())) { 
                    labIDs.put(l.getDBID(), l);
                    labNames.put(l.getName(), l);
                }
                
                return l;
            }            
            rs.close();
            return null;
        }
    }
    
    public Lab loadLab(int dbid) throws SQLException { 
        if(labIDs.containsKey(dbid)) { return labIDs.get(dbid); }

        Lab l = null;
        synchronized(loadLabs) {
            loadLabs.setInt(1, dbid);
            ResultSet rs = loadLabs.executeQuery();
            if(rs.next()) { 
                l = new Lab(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Lab DBID: " + dbid);
            }
        }        
        labIDs.put(dbid, l);
        labNames.put(l.getName(), l);
        return l;
    }
    
    public Collection<Lab> loadAllLabs(Collection<Integer> dbids) throws SQLException {

        LinkedList<Lab> values = new LinkedList<Lab>();
        for(int dbid : dbids) { values.addLast(loadLab(dbid)); }
        return values;
    }

    public Collection<Lab> loadAllLabs() throws SQLException {
        
        HashSet<Lab> values = new HashSet<Lab>();
        ResultSet rs = loadAllLabs.executeQuery();

        while(rs.next()) { 
        	Lab l = new Lab(rs);
            values.add(l);
            labNames.put(l.getName(), l);
            labIDs.put(l.getDBID(),l);
        }
        rs.close();
        return values;
    }	
	
    private int insertLab(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }

    
    //////////////////
    // CellLine stuff
    //////////////////

    public CellLine getCellLine(String name) throws SQLException { 
        synchronized (loadCellsByName) {
            loadCellsByName.setString(1, name);
            ResultSet rs = loadCellsByName.executeQuery();
            
            if(rs.next()) { 
                CellLine c = new CellLine(rs);
                rs.close();
                
                if(!cellIDs.containsKey(c.getDBID())) { 
                    cellIDs.put(c.getDBID(), c);
                    cellNames.put(c.getName(), c);
                }
                
                return c;
            }
            rs.close();
        }
        int id = insertCellLine(name);
        return loadCellLine(id);
    }
    
    public CellLine findCellLine(String name) throws SQLException { 
        synchronized (loadCellsByName) {
            loadCellsByName.setString(1, name);
            ResultSet rs = loadCellsByName.executeQuery();
            
            if(rs.next()) { 
                CellLine c = new CellLine(rs);
                rs.close();
                
                if(!cellIDs.containsKey(c.getDBID())) { 
                    cellIDs.put(c.getDBID(), c);
                    cellNames.put(c.getName(), c);
            }
                
                return c;
            }            
            rs.close();
            return null;
        }
    }
    
    public CellLine loadCellLine(int dbid) throws SQLException { 
        if(cellIDs.containsKey(dbid)) { return cellIDs.get(dbid); }

        CellLine c = null;
        synchronized(loadCells) {
            loadCells.setInt(1, dbid);
            ResultSet rs = loadCells.executeQuery();
            if(rs.next()) { 
                c = new CellLine(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Cells DBID: " + dbid);
            }
        }        
        cellIDs.put(dbid, c);
        cellNames.put(c.getName(), c);
        return c;
    }
    
    public Collection<CellLine> loadAllCellLines(Collection<Integer> dbids) throws SQLException {

        LinkedList<CellLine> values = new LinkedList<CellLine>();
        for(int dbid : dbids) { values.addLast(loadCellLine(dbid)); }
        return values;
    }

    public Collection<CellLine> loadAllCellLines() throws SQLException {
        
        HashSet<CellLine> values = new HashSet<CellLine>();
        ResultSet rs = loadAllCells.executeQuery();

        while(rs.next()) { 
            CellLine c = new CellLine(rs);
            values.add(c);
            cellNames.put(c.getName(), c);
            cellIDs.put(c.getDBID(),c);
        }
        rs.close();
        return values;
    }	

    private int insertCellLine(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }

    
    //////////////////
    // ExptCondition stuff
    //////////////////
    
    public ExptCondition getExptCondition(String name) throws SQLException { 
        synchronized (loadCondByName) {
            loadCondByName.setString(1, name);
            ResultSet rs = loadCondByName.executeQuery();
            
            if(rs.next()) { 
                ExptCondition c = new ExptCondition(rs);
                rs.close();
                
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
                
                return c;
            }
            rs.close();
        }
        int id = insertExptCondition(name);
        return loadExptCondition(id);
    }        

    public ExptCondition findExptCondition(String name) throws SQLException { 
        synchronized (loadCondByName) {
            loadCondByName.setString(1, name);
            ResultSet rs = loadCondByName.executeQuery();
            
            if(rs.next()) { 
                ExptCondition c = new ExptCondition(rs);
                rs.close();
                
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
                
                return c;
            }            
            rs.close();
            return null;
        }
    }
    
    public ExptCondition loadExptCondition(int dbid) throws SQLException { 
        if(condIDs.containsKey(dbid)) {  return condIDs.get(dbid); }
        
        ExptCondition c = null;
        synchronized(loadCond) {
            loadCond.setInt(1, dbid);
            ResultSet rs = loadCond.executeQuery();
            
            if(rs.next()) { 
                c = new ExptCondition(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ExptCondition DBID: " + dbid);
            }
        }
        
        condIDs.put(dbid, c);
        condNames.put(c.getName(), c);
        return c;
    }

    public Collection<ExptCondition> loadAllExptConditions(Collection<Integer> dbids) throws SQLException {

        LinkedList<ExptCondition> values = new LinkedList<ExptCondition>();
        for(int dbid : dbids) { values.addLast(loadExptCondition(dbid)); }
        return values;
    }

    public Collection<ExptCondition> loadAllExptConditions() throws SQLException { 
        HashSet<ExptCondition> values = new HashSet<ExptCondition>();
        ResultSet rs = loadAllCond.executeQuery();

        while(rs.next()) { 
            ExptCondition c = new ExptCondition(rs);
            values.add(c);
            condNames.put(c.getName(), c);
            condIDs.put(c.getDBID(),c);
        }
        rs.close();

        return values;
    }	

    private int insertExptCondition(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }

    
    //////////////////
    // ExptTarget stuff
    //////////////////
    
    public ExptTarget findExptTarget(String name) throws SQLException { 
        synchronized (loadTargetsByName) {
            loadTargetsByName.setString(1, name);
            ResultSet rs = loadTargetsByName.executeQuery();
            
            if(rs.next()) { 
                ExptTarget c = new ExptTarget(rs);
                rs.close();
                
                if(!targetIDs.containsKey(c.getDBID())) { 
                    targetIDs.put(c.getDBID(), c);
                    targetNames.put(c.getName(), c);
                }
                
                return c;
            }        
            rs.close();
        }
        return null;
    }
    
    public ExptTarget getExptTarget(String name) throws SQLException { 
        synchronized (loadTargetsByName) {
            loadTargetsByName.setString(1, name);
            ResultSet rs = loadTargetsByName.executeQuery();
            
            if(rs.next()) { 
                ExptTarget c = new ExptTarget(rs);
                rs.close();
                
                if(!targetIDs.containsKey(c.getDBID())) { 
                    targetIDs.put(c.getDBID(), c);
                    targetNames.put(c.getName(), c);
                }
                
                return c;
            } else { 
                rs.close();
            }
        }
        int id = insertFactor(name);
        return loadExptTarget(id);
    }
    
    public ExptTarget loadExptTarget(int dbid) throws SQLException { 
        if(targetIDs.containsKey(dbid)) { return targetIDs.get(dbid); }
        
        
        ExptTarget c = null;
        synchronized(loadTargets) {
            loadTargets.setInt(1, dbid);
            ResultSet rs = loadTargets.executeQuery();
            
            if(rs.next()) { 
                c = new ExptTarget(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ExptTarget DBID: " + dbid);
            }
        }
        
        targetIDs.put(dbid, c);
        targetNames.put(c.getName(), c);
        return c;
    }
  	
    public Collection<ExptTarget> loadAllExptTargets(Collection<Integer> dbids) throws SQLException {

        LinkedList<ExptTarget> values = new LinkedList<ExptTarget>();
        for(int dbid : dbids) { values.addLast(loadExptTarget(dbid)); }
        return values;
    }

    public Collection<ExptTarget> loadAllExptTargets() throws SQLException { 
        HashSet<ExptTarget> values = new HashSet<ExptTarget>();
        ResultSet rs = loadAllTargets.executeQuery();

        while(rs.next()) { 
            ExptTarget f = new ExptTarget(rs);
            values.add(f);
            targetNames.put(f.getName(), f);
            targetIDs.put(f.getDBID(),f);
        }
        rs.close();
        return values;
    }	
	
    private int insertFactor(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }

    //////////////////
    // ExptType stuff
    //////////////////
    
    public ExptType getExptType(String name) throws SQLException { 
        synchronized (loadExptTypesByName) {
            loadExptTypesByName.setString(1, name);
            ResultSet rs = loadExptTypesByName.executeQuery();
            
            if(rs.next()) { 
            	ExptType e = new ExptType(rs);
                rs.close();
                
                if(!exptTypeIDs.containsKey(e.getDBID())) { 
                    exptTypeIDs.put(e.getDBID(), e);
                    exptTypeNames.put(e.getName(), e);
                }
                
                return e;
            }
            rs.close();
        }
        int id = insertExptType(name);
        return loadExptType(id);
    }
    
    public ExptType findExptType(String name) throws SQLException { 
        synchronized (loadExptTypesByName) {
            loadExptTypesByName.setString(1, name);
            ResultSet rs = loadExptTypesByName.executeQuery();
            
            if(rs.next()) { 
            	ExptType e = new ExptType(rs);
                rs.close();
                
                if(!exptTypeIDs.containsKey(e.getDBID())) { 
                    exptTypeIDs.put(e.getDBID(), e);
                    exptTypeNames.put(e.getName(), e);
                }
                return e;
            }            
            rs.close();
            return null;
        }
    }
    
    public ExptType loadExptType(int dbid) throws SQLException { 
        if(exptTypeIDs.containsKey(dbid)) { return exptTypeIDs.get(dbid); }

        ExptType e = null;
        synchronized(loadExptTypes) {
            loadExptTypes.setInt(1, dbid);
            ResultSet rs = loadExptTypes.executeQuery();
            if(rs.next()) { 
                e = new ExptType(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ExptType DBID: " + dbid);
            }
        }        
        exptTypeIDs.put(dbid, e);
        exptTypeNames.put(e.getName(), e);
        return e;
    }
    
    public Collection<ExptType> loadAllExptTypes(Collection<Integer> dbids) throws SQLException {

        LinkedList<ExptType> values = new LinkedList<ExptType>();
        for(int dbid : dbids) { values.addLast(loadExptType(dbid)); }
        return values;
    }

    public Collection<ExptType> loadAllExptTypes() throws SQLException {
        
        HashSet<ExptType> values = new HashSet<ExptType>();
        ResultSet rs = loadAllExptTypes.executeQuery();

        while(rs.next()) { 
        	ExptType e = new ExptType(rs);
            values.add(e);
            exptTypeNames.put(e.getName(), e);
            exptTypeIDs.put(e.getDBID(),e);
        }
        rs.close();
        return values;
    }	
	
    private int insertExptType(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }

    //////////////////
    // ReadType stuff
    //////////////////
    
    public ReadType getReadType(String name) throws SQLException { 
        synchronized (loadReadTypesByName) {
            loadReadTypesByName.setString(1, name);
            ResultSet rs = loadReadTypesByName.executeQuery();
            
            if(rs.next()) { 
            	ReadType r = new ReadType(rs);
                rs.close();
                
                if(!readTypeIDs.containsKey(r.getDBID())) { 
                    readTypeIDs.put(r.getDBID(), r);
                    readTypeNames.put(r.getName(), r);
                }
                return r;
            }
            rs.close();
        }
        int id = insertReadType(name);
        return loadReadType(id);
    }
    
    public ReadType findReadType(String name) throws SQLException { 
        synchronized (loadReadTypesByName) {
            loadReadTypesByName.setString(1, name);
            ResultSet rs = loadReadTypesByName.executeQuery();
            
            if(rs.next()) { 
            	ReadType r = new ReadType(rs);
                rs.close();
                
                if(!readTypeIDs.containsKey(r.getDBID())) { 
                    readTypeIDs.put(r.getDBID(), r);
                    readTypeNames.put(r.getName(), r);
                }
                return r;
            }            
            rs.close();
            return null;
        }
    }
    
    public ReadType loadReadType(int dbid) throws SQLException { 
        if(readTypeIDs.containsKey(dbid)) { return readTypeIDs.get(dbid); }

        ReadType e = null;
        synchronized(loadReadTypes) {
            loadReadTypes.setInt(1, dbid);
            ResultSet rs = loadReadTypes.executeQuery();
            if(rs.next()) { 
                e = new ReadType(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown ReadType DBID: " + dbid);
            }
        }        
        readTypeIDs.put(dbid, e);
        readTypeNames.put(e.getName(), e);
        return e;
    }
    
    public Collection<ReadType> loadAllReadTypes(Collection<Integer> dbids) throws SQLException {

        LinkedList<ReadType> values = new LinkedList<ReadType>();
        for(int dbid : dbids) { values.addLast(loadReadType(dbid)); }
        return values;
    }

    public Collection<ReadType> loadAllReadTypes() throws SQLException {
        
        HashSet<ReadType> values = new HashSet<ReadType>();
        ResultSet rs = loadAllReadTypes.executeQuery();

        while(rs.next()) { 
        	ReadType r = new ReadType(rs);
            values.add(r);
            readTypeNames.put(r.getName(), r);
            readTypeIDs.put(r.getDBID(),r);
        }
        rs.close();
        return values;
    }	
	
    private int insertReadType(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }


    //////////////////
    // AlignType stuff
    //////////////////
    
    public AlignType getAlignType(String name) throws SQLException { 
        synchronized (loadAlignTypesByName) {
            loadAlignTypesByName.setString(1, name);
            ResultSet rs = loadAlignTypesByName.executeQuery();
            
            if(rs.next()) { 
            	AlignType a = new AlignType(rs);
                rs.close();
                
                if(!alignTypeIDs.containsKey(a.getDBID())) { 
                    alignTypeIDs.put(a.getDBID(), a);
                    alignTypeNames.put(a.getName(), a);
                }
                return a;
            }
            rs.close();
        }
        int id = insertAlignType(name);
        return loadAlignType(id);
    }
    
    public AlignType findAlignType(String name) throws SQLException { 
        synchronized (loadAlignTypesByName) {
            loadAlignTypesByName.setString(1, name);
            ResultSet rs = loadAlignTypesByName.executeQuery();
            
            if(rs.next()) { 
            	AlignType a = new AlignType(rs);
                rs.close();
                
                if(!alignTypeIDs.containsKey(a.getDBID())) { 
                    alignTypeIDs.put(a.getDBID(), a);
                    alignTypeNames.put(a.getName(), a);
                }
                return a;
            }            
            rs.close();
            return null;
        }
    }
    
    public AlignType loadAlignType(int dbid) throws SQLException { 
        if(alignTypeIDs.containsKey(dbid)) { return alignTypeIDs.get(dbid); }

        AlignType a = null;
        synchronized(loadAlignTypes) {
            loadAlignTypes.setInt(1, dbid);
            ResultSet rs = loadAlignTypes.executeQuery();
            if(rs.next()) { 
                a = new AlignType(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown AlignType DBID: " + dbid);
            }
        }        
        alignTypeIDs.put(dbid, a);
        alignTypeNames.put(a.getName(), a);
        return a;
    }
    
    public Collection<AlignType> loadAllAlignTypes(Collection<Integer> dbids) throws SQLException {

        LinkedList<AlignType> values = new LinkedList<AlignType>();
        for(int dbid : dbids) { values.addLast(loadAlignType(dbid)); }
        return values;
    }

    public Collection<AlignType> loadAllAlignTypes() throws SQLException {
        
        HashSet<AlignType> values = new HashSet<AlignType>();
        ResultSet rs = loadAllReadTypes.executeQuery();

        while(rs.next()) { 
        	AlignType a = new AlignType(rs);
            values.add(a);
            alignTypeNames.put(a.getName(), a);
            alignTypeIDs.put(a.getDBID(),a);
        }
        rs.close();
        return values;
    }	
	
    private int insertAlignType(String n) throws SQLException {
    	Statement s = null;
    	ResultSet rs = null;
    	int id=-1;
    	try{
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
	    } finally {
	        if (rs != null) {
	            try {
	                rs.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }

	        if (s != null) {
	            try {
	                s.close();
	            } catch (SQLException ex) {
	                // ignore
	            }
	        }
	    }
        return id;
    }


}
