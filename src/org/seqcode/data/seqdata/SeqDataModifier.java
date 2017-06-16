package org.seqcode.data.seqdata;

import java.io.IOException;
import java.security.AccessControlException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashSet;
import java.util.Set;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.core.CellLine;
import org.seqcode.data.core.ExptCondition;
import org.seqcode.data.core.ExptTarget;
import org.seqcode.data.core.Lab;
import org.seqcode.data.core.MetadataLoader;
import org.seqcode.data.core.MetadataModifier;
import org.seqcode.data.readdb.ACLChangeEntry;
import org.seqcode.data.readdb.Client;
import org.seqcode.data.readdb.ClientException;
import org.seqcode.gseutils.NotFoundException;

/**
 * SeqDataModifier collects interactions that modify the seqdata database
 * 
 * @author mahony
 *
 */
public class SeqDataModifier {
	public static String role = "seqdata";

	private SeqDataLoader seqLoader;
	private Client client = null;
	private MetadataModifier metaModifier = null;
	private boolean closed = true;

	public SeqDataModifier(SeqDataLoader loader) throws AccessControlException, SQLException {
		seqLoader = loader;
		client = seqLoader.getClient();
		metaModifier = new MetadataModifier();
		if (!seqLoader.getMyUser().isAdmin()) {
			throw new AccessControlException("SeqDataModifier: only admins can modify seqdata!");
		}
		closed = false;
	}

	public void deleteAlignmentParameters(SeqAlignment align) throws SQLException {
		Connection cxn = null;
		Statement del = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			del = cxn.createStatement();
			del.execute("delete from alignmentparameters where alignment = " + align.getDBID());
			cxn.commit();
		} catch (SQLException e) {
			cxn.rollback();
			throw new DatabaseException(e.toString(), e);
		} finally {
			if (del != null) {
				try {
					del.close();
				} catch (SQLException ex) {
				}
			}
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void deleteSeqAlignment(SeqAlignment align) throws SQLException {
		Connection cxn = null;
		PreparedStatement ps = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			ps = SeqAlignment.createDeleteByIDStatement(cxn);
			ps.setInt(1, align.getDBID());
			ps.execute();
			cxn.commit();
		} catch (SQLException e) {
			cxn.rollback();
			throw new DatabaseException(e.toString(), e);
		} finally {
			if (ps != null) {
				try {
					ps.close();
				} catch (SQLException ex) {
				}
			}
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void deleteSeqExpt(SeqExpt expt) throws SQLException {
		Connection cxn = null;
		PreparedStatement ps = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			Lab lab = expt.getLab();
			ExptCondition cond = expt.getExptCondition();
			ExptTarget target = expt.getExptTarget();
			CellLine cells = expt.getCellLine();

			ps = SeqExpt.createDeleteByDBID(cxn);
			ps.setInt(1, expt.getDBID());
			ps.execute();

			// Delete core.lab if no SeqExpts depend
			if (seqLoader.loadExperiments(lab).size() == 0)
				deleteLab(lab);
			// Delete core.exptcondition if no SeqExpts depend
			if (seqLoader.loadExperiments(cond).size() == 0)
				deleteExptCondition(cond);
			// Delete core.expttarget if no SeqExpts depend
			if (seqLoader.loadExperiments(target).size() == 0)
				deleteExptTarget(target);
			// Delete core.cellline if no SeqExpts depend
			if (seqLoader.loadExperiments(cells).size() == 0)
				deleteCellLine(cells);
			cxn.commit();
		} catch (SQLException e) {
			cxn.rollback();
			throw new DatabaseException(e.toString(), e);
		} finally {
			if (ps != null) {
				try {
					ps.close();
				} catch (SQLException ex) {
				}
			}
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void deleteLab(Lab lab) throws SQLException {
		System.err.println("Deleting lab: " + lab.getName() + "\t" + lab.getDBID());
		metaModifier.deleteLab(lab.getDBID());
	}

	public void deleteExptCondition(ExptCondition cond) throws SQLException {
		System.err.println("Deleting condition: " + cond.getName() + "\t" + cond.getDBID());
		metaModifier.deleteCond(cond.getDBID());
	}

	public void deleteExptTarget(ExptTarget target) throws SQLException {
		System.err.println("Deleting target: " + target.getName() + "\t" + target.getDBID());
		metaModifier.deleteTarget(target.getDBID());
	}

	public void deleteCellLine(CellLine cells) throws SQLException {
		System.err.println("Deleting cell-line: " + cells.getName() + "\t" + cells.getDBID());
		metaModifier.deleteCell(cells.getDBID());
	}

	public void coreCleanup() throws SQLException {
		for (Lab lab : seqLoader.getMetadataLoader().loadAllLabs(true)) {
			// Delete core.lab if no SeqExpts depend
			if (seqLoader.loadExperiments(lab).size() == 0)
				deleteLab(lab);
		}
		for (ExptCondition cond : seqLoader.getMetadataLoader().loadAllExptConditions(true)) {
			// Delete core.exptcondition if no SeqExpts depend
			if (seqLoader.loadExperiments(cond).size() == 0)
				deleteExptCondition(cond);
		}
		for (ExptTarget target : seqLoader.getMetadataLoader().loadAllExptTargets(true)) {
			// Delete core.expttarget if no SeqExpts depend
			if (seqLoader.loadExperiments(target).size() == 0)
				deleteExptTarget(target);
		}
		for (CellLine cells : seqLoader.getMetadataLoader().loadAllCellLines(true)) {
			// Delete core.cellline if no SeqExpts depend
			if (seqLoader.loadExperiments(cells).size() == 0)
				deleteCellLine(cells);
		}
	}

	public void updateSeqExpt(MetadataLoader mloader, SeqExpt expt, String updateExptType, String updateLab,
			String updateCond, String updateTarget, String updateCell, String updateRep, String updatePubSrc,
			String updatePubID, String updateCollabExptID) throws SQLException, DuplicateDatabaseEntryException {
		Connection cxn = null;
		PreparedStatement update = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(true);
			String updateName = updateLab + " " + updateCond + " " + updateTarget + " " + updateCell;

			SeqExpt testExpt = seqLoader.findExperiment(updateName, updateRep);
			if (testExpt != null && testExpt.getDBID() != expt.getDBID()) // It's
																			// okay
																			// if
																			// these
																			// are
																			// the
																			// same
																			// experiments
																			// (you
																			// might
																			// sometimes
																			// want
																			// to
																			// just
																			// update
																			// the
																			// publication
																			// source,
																			// etc).
				throw new DuplicateDatabaseEntryException(
						"SeqDataModifier.updateSeqExpt wants to create a duplicate SeqExpt");
			else {
				update = SeqExpt.createShortUpdateWithID(cxn);
				update.setString(1, updateName);
				update.setString(2, updateRep);
				update.setInt(3, expt.getOrganism().getDBID());
				update.setInt(4, mloader.loadExptType(updateExptType, true, false).getDBID());
				update.setInt(5, mloader.loadLab(updateLab, true, false).getDBID());
				update.setInt(6, mloader.loadExptCondition(updateCond, true, false).getDBID());
				update.setInt(7, mloader.loadExptTarget(updateTarget, true, false).getDBID());
				update.setInt(8, mloader.loadCellLine(updateCell, true, false).getDBID());
				update.setString(9, updateCollabExptID);
				update.setString(10, updatePubSrc);
				update.setString(11, updatePubID);
				update.setInt(12, expt.getDBID());
				update.execute();
				update.close();

				try {
					SeqExpt testExpt2 = seqLoader.loadExperiment(updateName, updateRep);
				} catch (NotFoundException e2) {
					// failed again means the insert failed. you lose
					throw new DatabaseException("Couldn't update experiment for " + updateName + "," + updateRep);
				}
			}
		} catch (SQLException e) {
			throw new DatabaseException(e.toString(), e);
		} finally {
			if (update != null) {
				try {
					update.close();
				} catch (SQLException ex) {
				}
			}
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void updateSeqAlignmentHitCounts(SeqAlignment align, Integer singlecount, Float singleweight,
			Integer singletype2count, Float singletype2weight, Integer paircount, Float pairweight)
			throws SQLException {
		Connection cxn = null;
		PreparedStatement update = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			int id = align.getDBID();
			update = SeqAlignment.createUpdateHitsAndWeights(cxn);
			System.err.println("Updating counts for alignment: " + id + " (" + align.getName() + ")");
			System.err.println("\tnumhits=" + singlecount);
			System.err.println("\ttotalweight=" + singleweight);
			System.err.println("\tnumtype2hits=" + singletype2count);
			System.err.println("\ttotaltype2weight=" + singletype2weight);
			System.err.println("\tnumpairs=" + paircount);
			System.err.println("\ttotalpairweight=" + pairweight);
			update.setInt(1, singlecount);
			update.setFloat(2, singleweight);
			update.setInt(3, singletype2count);
			update.setFloat(4, singletype2weight);
			update.setInt(5, paircount);
			update.setFloat(6, pairweight);
			update.setInt(7, id);
			update.execute();
			cxn.commit();
		} catch (SQLException e) {
			cxn.rollback();
			throw new DatabaseException(e.toString(), e);
		} finally {
			if (update != null) {
				try {
					update.close();
				} catch (SQLException ex) {
				}
			}
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	public void updateSeqAlignmentPermissions(SeqAlignment align, String permissions) throws SQLException {
		Connection cxn = null;
		PreparedStatement permUpdate = null;
		try {
			cxn = DatabaseConnectionManager.getConnection(role);
			cxn.setAutoCommit(false);
			permUpdate = SeqAlignment.createUpdatePermissions(cxn);
			permUpdate.setString(1, permissions);
			permUpdate.setInt(2, align.getDBID());
			permUpdate.execute();
			permUpdate.close();
			cxn.commit();
		} catch (SQLException e) {
			cxn.rollback();
			throw new DatabaseException(e.toString(), e);
		} finally {
			if (permUpdate != null) {
				try {
					permUpdate.close();
				} catch (SQLException ex) {
				}
			}
			if (cxn != null)
				try {
					cxn.close();
				} catch (Exception ex) {
					throw new DatabaseException("Couldn't close connection with role " + role, ex);
				}
		}
	}

	/**
	 * Update the permissions for a SeqAlignment SeqAlignment align: alignment
	 * to change String princ : user name String op : operation [add|delete]
	 * String acl [read|write|admin]
	 * 
	 * @param princ
	 */
	public void changeAlignmentACL(SeqAlignment align, String princ, String op, String acl) {
		Set<ACLChangeEntry> changes = new HashSet<ACLChangeEntry>();
		changes.add(new ACLChangeEntry(ACLChangeEntry.opCode(op), ACLChangeEntry.aclCode(acl), princ));
		try {
			client.setACL(new Integer(align.getDBID()).toString(), changes);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClientException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Update multiple permissions for a SeqAlignment SeqAlignment align:
	 * alignment to change String[] princs : user name String[] ops : operation
	 * [add|delete] String[] acls [read|write|admin]
	 * 
	 * @param princ
	 */
	public void changeAlignmentACLmulti(SeqAlignment align, String[] princs, String[] ops, String[] acls) {
		if (princs.length == ops.length && ops.length == acls.length) {
			Set<ACLChangeEntry> changes = new HashSet<ACLChangeEntry>();
			for (int i = 0; i < princs.length; i++)
				changes.add(
						new ACLChangeEntry(ACLChangeEntry.opCode(ops[i]), ACLChangeEntry.aclCode(acls[i]), princs[i]));
			try {
				client.setACL(new Integer(align.getDBID()).toString(), changes);
			} catch (IOException e) {
				e.printStackTrace();
			} catch (ClientException e) {
				e.printStackTrace();
			}
		} else {
			System.err.println("changeAlignmentACLmulti: input arrays should be the same lengths");
		}
	}

	public class DuplicateDatabaseEntryException extends Exception {
		public DuplicateDatabaseEntryException(String message) {
			super(message);
		}
	}

}
