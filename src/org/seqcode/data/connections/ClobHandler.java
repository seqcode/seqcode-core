package org.seqcode.data.connections;

import oracle.jdbc.OraclePreparedStatement;
import java.sql.*;

public class ClobHandler {

	public static void setClob(Connection cxn, PreparedStatement ps, int index, String clobValue) throws SQLException {
		if (DatabaseConnectionManager.isOracle(cxn)) {
			OraclePreparedStatement ops = (OraclePreparedStatement) ps;
			ops.setStringForClob(index, clobValue);
		} else {
			ps.setString(index, clobValue);
		}
	}
}
