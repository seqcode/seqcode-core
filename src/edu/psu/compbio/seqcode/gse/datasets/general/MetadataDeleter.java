package edu.psu.compbio.seqcode.gse.datasets.general;

import java.sql.PreparedStatement;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;


public class MetadataDeleter implements edu.psu.compbio.seqcode.gse.utils.Closeable {
    
    public static final String role = "core";
	
    private java.sql.Connection cxn;
    
    private PreparedStatement deleteLabs, deleteCells, deleteCond, deleteTargets;
	
    public MetadataDeleter() throws SQLException { 
        try {
            cxn = DatabaseFactory.getConnection(role);
            cxn.setAutoCommit(false);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }

        deleteLabs = cxn.prepareStatement("delete from lab where id=?");
        deleteCells = cxn.prepareStatement("delete from cellline where id=?");
        deleteCond = cxn.prepareStatement("delete from exptcondition where id=?");
        deleteTargets = cxn.prepareStatement("delete from expttarget where id=?");
        
    }
	
    public boolean isClosed() { return cxn==null; }

    public void close() { 
        try {
        	deleteLabs.close(); deleteLabs=null;
        	deleteCells.close(); deleteCells = null;
        	deleteCond.close();  deleteCond = null;
        	deleteTargets.close(); deleteTargets = null;
        } catch (SQLException e) {
            e.printStackTrace();
        }
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }
    
    public java.sql.Connection getConnection() { return cxn; }
 
   
    public void deleteLab(int dbid) throws SQLException { 
    	synchronized(deleteLabs) {
    		deleteLabs.setInt(1, dbid);
            deleteLabs.execute();
            deleteLabs.close();
            cxn.commit();
        }        
    }
    
    public void deleteCell(int dbid) throws SQLException { 
    	synchronized(deleteCells) {
    		deleteCells.setInt(1, dbid);
            deleteCells.execute();
            deleteCells.close();
            cxn.commit();
        }        
    }
    
    public void deleteCond(int dbid) throws SQLException { 
    	synchronized(deleteCond) {
    		deleteCond.setInt(1, dbid);
            deleteCond.execute();
            deleteCond.close();
            cxn.commit();
        }        
    }
    
    public void deleteTarget(int dbid) throws SQLException { 
    	synchronized(deleteTargets) {
    		deleteTargets.setInt(1, dbid);
            deleteTargets.execute();
            deleteTargets.close();
            cxn.commit();
        }        
    }
 
}
