/**
 */
package org.seqcode.gsebricks.verbs;

import java.util.*;

import org.seqcode.gsebricks.iterators.EmptyIterator;

import java.sql.*;


/**
 * @author tdanford
 */
public class TableXRef implements Expander<String,String>, org.seqcode.gseutils.Closeable {
	
	private Connection cxn;
	private PreparedStatement ps;
	private String tableName, matchField, returnField;
	
	public TableXRef(Connection c, String t, String match, String ret) 
		throws SQLException {
		
		cxn = c;
		tableName = t;
		matchField = match;
		returnField = ret;
		
		ps = cxn.prepareStatement("select " + returnField + " from " + tableName + " " +
				"where " + matchField + "=?");
	}

	public Iterator<String> execute(String input) {
		try {
			ps.setString(1, input);
			ResultSet rs = ps.executeQuery();
			LinkedList<String> results = new LinkedList<String>();
			while(rs.next()) { 
				results.add(rs.getString(1));
			}
			rs.close();
			return results.iterator();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return new EmptyIterator<String>();
	}

	public void close() {
		try {
			ps.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		ps = null;
		cxn = null;
	}

	public boolean isClosed() {
		return cxn==null;
	}
}
