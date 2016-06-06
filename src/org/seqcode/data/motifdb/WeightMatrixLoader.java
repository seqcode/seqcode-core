package org.seqcode.data.motifdb;

import java.sql.*;
import java.util.*;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.UnknownRoleException;
import org.seqcode.genome.Species;
import org.seqcode.utils.NotFoundException;


public class WeightMatrixLoader implements org.seqcode.utils.Closeable {

    public WeightMatrixLoader() {}

    private Collection<String> queryWMTable(String field) {
        try {
            java.sql.Connection cxn = 
                DatabaseConnectionManager.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select unique(" + field + ") from weightmatrix order by " + field);
            ResultSet rs = ps.executeQuery();
            ArrayList<String> out = new ArrayList<String>();
            while (rs.next()) {
                if (rs.getString(1) == null) {
                    System.err.println("Found a null string in " + field);
                } else {
                    out.add(rs.getString(1));
                }
            }
            rs.close();
            ps.close();
            if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role annotations", ex); }
            return out;
        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        }
    }
    public Collection<String> getNames() {
        return queryWMTable("name");
    }
    public Collection<String> getVersions() {
        return queryWMTable("version");
    }

    public Collection<String> getTypes() {
        return queryWMTable("type");
    }

    public WeightMatrix query(int speciesid,
                              String name,
                              String version) {
        WeightMatrix wm = null;
        try {            
            java.sql.Connection cxn = 
                DatabaseConnectionManager.getConnection("annotations");
            String query = "select wm.id from weightmatrix wm where wm.species = ? and " +
                    " wm.name = ? and wm.version = ?";
            PreparedStatement ps = cxn.prepareStatement(query);
            ps.setInt(1,speciesid);
            ps.setString(2,name);
            ps.setString(3,version);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                wm = WeightMatrix.getWeightMatrix(rs.getInt(1));
            }
            rs.close();
            ps.close();
            if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role annotations", ex); }
            if (wm != null) {
                return wm;
            } else {
                throw new NotFoundException ("Couldn't find " + name + "," + version + " in species " + speciesid);
            }

        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        }
    
    }
    
    public Collection<WeightMatrix> loadMatrices(Species species) throws SQLException { 
    	int speciesID = species.getDBID();
        String query = "select m.id, m.species, m.name, m.version, m.type, m.bg_model_map_id, c.position, c.letter, c.weight from weightmatrix m, " +
            " weightmatrixcols c where m.id = c.weightmatrix and m.species = ? order by c.weightmatrix, c.position desc";

        java.sql.Connection cxn = 
            DatabaseConnectionManager.getConnection("annotations");
    	PreparedStatement wmStatement = cxn.prepareStatement(query);
    	
    	wmStatement.setInt(1, speciesID);
    	ResultSet wmResults = wmStatement.executeQuery();
        Collection<WeightMatrix> matrices = WeightMatrix.getWeightMatrices(wmResults);
    	wmResults.close();
    	wmStatement.close();
    	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role annotations", ex); }
    	
    	return matrices;
    }

    public Collection<WeightMatrix> query(String name,
                                          String version,
                                          String type) {
        if (name == null && version == null && type == null) {
            return WeightMatrix.getAllWeightMatrices();
        }
        try {
            java.sql.Connection cxn = 
                DatabaseConnectionManager.getConnection("annotations");
            ArrayList<WeightMatrix> out = new ArrayList<WeightMatrix>();
            String query;
            query = "select wm.id, wm.species, wm.name, wm.version, wm.type, wm.bg_model_map_id, c.position, c.letter, c.weight " +
                " from weightmatrix wm, weightmatrixcols c where wm.id = c.weightmatrix ";
            if (name != null) {
                query += " and wm.name = ? ";
            }
            if (version != null) {
                query += " and wm.version = ? ";
                }
            if (type != null) {
                query += " and wm.type = ? ";
            }
            query += " order by c.weightmatrix, c.position desc";
            PreparedStatement ps = cxn.prepareStatement(query);
            int index = 1;
            if (name != null) {
                ps.setString(index++,name);
            }
            if (version != null) {
                ps.setString(index++,version);
            }
            if (type != null) {
                ps.setString(index++,type);
            }
            ResultSet rs = ps.executeQuery();
            Collection<WeightMatrix> matrices = WeightMatrix.getWeightMatrices(rs);
            rs.close();
            ps.close();
            if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role annotations", ex); }

            return matrices;
        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        }
    }

    public void close() {}
    public boolean isClosed() {return true;}
}
