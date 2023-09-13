package org.seqcode.data.connections;

import java.sql.*;

public class ClobHandler {

	public static void setClob(Connection cxn, PreparedStatement ps, int index, String clobValue) throws SQLException { 
		if(DatabaseConnectionManager.isOracle(cxn)) { 
			
		} else { 
			ps.setString(index, clobValue);
		}
	}
}
