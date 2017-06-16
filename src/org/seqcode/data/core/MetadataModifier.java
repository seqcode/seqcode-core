package org.seqcode.data.core;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.UnknownRoleException;

public class MetadataModifier {

	public static final String role = "core";

	public MetadataModifier() throws SQLException {
	}

	public void deleteLab(int dbid) throws SQLException {
		Connection cxn = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			PreparedStatement deleteLabs = cxn.prepareStatement("delete from lab where id=?");
			deleteLabs.setInt(1, dbid);
			deleteLabs.execute();
			deleteLabs.close();
			cxn.commit();
		} catch (UnknownRoleException e) {
			cxn.rollback();
			throw new IllegalArgumentException("Unknown role: " + role, e);
		} finally {
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void deleteCell(int dbid) throws SQLException {
		Connection cxn = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			PreparedStatement deleteCells = cxn.prepareStatement("delete from cellline where id=?");
			deleteCells.setInt(1, dbid);
			deleteCells.execute();
			deleteCells.close();
			cxn.commit();
		} catch (UnknownRoleException e) {
			cxn.rollback();
			throw new IllegalArgumentException("Unknown role: " + role, e);
		} finally {
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void deleteCond(int dbid) throws SQLException {
		Connection cxn = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			PreparedStatement deleteCond = cxn.prepareStatement("delete from exptcondition where id=?");
			deleteCond.setInt(1, dbid);
			deleteCond.execute();
			deleteCond.close();
			cxn.commit();
		} catch (UnknownRoleException e) {
			cxn.rollback();
			throw new IllegalArgumentException("Unknown role: " + role, e);
		} finally {
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void deleteTarget(int dbid) throws SQLException {
		Connection cxn = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			PreparedStatement deleteTargets = cxn.prepareStatement("delete from expttarget where id=?");
			deleteTargets.setInt(1, dbid);
			deleteTargets.execute();
			deleteTargets.close();
			cxn.commit();
		} catch (UnknownRoleException e) {
			cxn.rollback();
			throw new IllegalArgumentException("Unknown role: " + role, e);
		} finally {
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

}
