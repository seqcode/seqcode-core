package org.seqcode.data.seqdata.tools;

import java.io.*;
import java.sql.*;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.UnknownRoleException;
import org.seqcode.data.core.*;
import org.seqcode.data.seqdata.*;
import org.seqcode.genome.Genome;
import org.seqcode.gseutils.*;

import com.mysql.jdbc.exceptions.jdbc4.MySQLIntegrityConstraintViolationException;

/**
 * Creates an experiment (if necessary) in the database and prints the DBID on
 * stdout. Use this with the readdb importer to load metadata to the mysql
 * database. I split this out from a combined experiment/alignment creator
 * because of what seems like simultaneous write conflicts during experiment
 * creation.
 *
 * Usage: CreateExpt --species "$SC;SGDv1" --align "name;replicate;alignment
 * version" --lab "Pugh" --expttype "CHIPSEQ" --expttarget "Gcn4" --cellline
 * "FY4" --exptcondition "YPD" --readtype "SINGLE" --aligntype "SINGLE"
 * --paramsfile params.txt --readlength 36
 */
public class CreateExpt {

	public static void main(String args[]) throws SQLException, IOException, NotFoundException {

		if (args.length == 0) {
			System.out.println("CreateExpt:\n" + "\t--species <species;genome>\n"
					+ "\t--align <name;replicate;version>\n" + "\t--expttype <CHIPSEQ/CHIPEXO/RNASEQ/etc>\n"
					+ "\t--lab <name>\n" + "\t--exptcondition <condition>\n" + "\t--expttarget <target>\n"
					+ "\t--cellline <cell line>\n" + "\t--readtype <SINGLE/PAIRED>\n" + "\t--readlength <int>\n"
					+ "\t--numreads <int>\n" + "\t--collabid <expt ID>\n" + "\t--publicsource <PMID/UNPUB>\n"
					+ "\t--publicdbid <GEO ID>\n" + "\t--fqfile <FQ filename>\n" + "\t--exptnote <notes about expt>\n");
		} else {

			java.sql.Connection cxn = null;
			PreparedStatement insert = null;
			PreparedStatement update = null;

			SeqExpt expt = null;
			SeqDataLoader loader = new SeqDataLoader();
			MetadataLoader core = loader.getMetadataLoader();
			boolean newExpt = false;

			try {
				cxn = DatabaseConnectionManager.getConnection("seqdata");
				cxn.setAutoCommit(true);
				Genome genome = Args.parseGenome(args).cdr();
				String alignname = Args.parseString(args, "align", null);
				String alignpieces[] = alignname.split(";");
				String etypestring = Args.parseString(args, "expttype", null);
				String labstring = Args.parseString(args, "lab", null);
				String conditionstring = Args.parseString(args, "exptcondition", null);
				String targetstring = Args.parseString(args, "expttarget", null);
				String cellsstring = Args.parseString(args, "cellline", null);
				String rtypestring = Args.parseString(args, "readtype", null);
				int readlength = Args.parseInteger(args, "readlength", 36);
				int numreads = Args.parseInteger(args, "numreads", 0);
				String collabid = Args.parseString(args, "collabid", null);
				String publicsource = Args.parseString(args, "publicsource", "UNPUB");
				String publicdbid = Args.parseString(args, "publicdbid", "NA");
				String fqfile = Args.parseString(args, "fqfile", null);
				String exptnote = Args.parseString(args, "exptnote", null);

				// SEQEXPERIMENT
				try {
					expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
				} catch (NotFoundException e) {
					try {
						// NotFound = create new experiment
						newExpt = true;
						System.err.println(
								"Creating experiment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
						insert = SeqExpt.createInsert(cxn);
						insert.setString(1, alignpieces[0]);
						insert.setString(2, alignpieces[1]);
						insert.setInt(3, genome.getSpeciesDBID());
						insert.setInt(4, core.loadExptType(etypestring, true, false).getDBID());
						insert.setInt(5, core.loadLab(labstring, true, false).getDBID());
						insert.setInt(6, core.loadExptCondition(conditionstring, true, false).getDBID());
						insert.setInt(7, core.loadExptTarget(targetstring, true, false).getDBID());
						insert.setInt(8, core.loadCellLine(cellsstring, true, false).getDBID());
						insert.setInt(9, core.loadReadType(rtypestring, true, false).getDBID());
						insert.setInt(10, readlength);
						insert.setInt(11, numreads);
						insert.setString(12, collabid);
						insert.setString(13, publicsource);
						insert.setString(14, publicdbid);
						insert.setString(15, fqfile);
						insert.setString(16, exptnote);
						insert.execute();
						insert.close();
					} catch (MySQLIntegrityConstraintViolationException ex) {
						// Duplicate primary keys - someone else added the same
						// experiment at the same time
						// Do nothing - it will try loading the experiment again
						// in the next step
					}
					try {
						expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
					} catch (NotFoundException e2) {
						/* failed again means the insert failed. you lose */
						// cxn.rollback(); //Can't use rollback with
						// autocommit=true?
						throw new DatabaseException("Couldn't create " + alignpieces[0] + "," + alignpieces[1]);
					}
				}
				if (!newExpt) {
					// Experiment exists: Update the old experiment
					// Check if the experiment information is the same before
					// bothering to update the database.
					if (!(expt.getName().equals(alignpieces[0]) && expt.getReplicate().equals(alignpieces[1])
							&& expt.getExptType().getDBID() == core.loadExptType(etypestring, true, false).getDBID()
							&& expt.getLab().getDBID() == core.loadLab(labstring, true, false).getDBID()
							&& expt.getExptCondition().getDBID() == core.loadExptCondition(conditionstring, true, false)
									.getDBID()
							&& expt.getExptTarget().getDBID() == core.loadExptTarget(targetstring, true, false)
									.getDBID()
							&& expt.getCellLine().getDBID() == core.loadCellLine(cellsstring, true, false).getDBID()
							&& expt.getReadType().getDBID() == core.loadReadType(rtypestring, true, false).getDBID()
							&& expt.getReadLength() == readlength && expt.getNumRead() == numreads
							&& expt.getCollabID().equals(collabid) && expt.getPublicSource().equals(publicsource)
							&& expt.getPublicDBID().equals(publicdbid) && expt.getFQFile().equals(fqfile)
							&& expt.getExptNote().equals(exptnote))) {

						int eID = expt.getDBID();
						System.err.println("Updating experiment " + eID + " " + alignpieces[0] + ";" + alignpieces[1]
								+ ";" + alignpieces[2]);
						update = SeqExpt.createUpdateWithID(cxn);
						update.setString(1, alignpieces[0]);
						update.setString(2, alignpieces[1]);
						update.setInt(3, genome.getSpeciesDBID());
						update.setInt(4, core.loadExptType(etypestring, true, false).getDBID());
						update.setInt(5, core.loadLab(labstring, true, false).getDBID());
						update.setInt(6, core.loadExptCondition(conditionstring, true, false).getDBID());
						update.setInt(7, core.loadExptTarget(targetstring, true, false).getDBID());
						update.setInt(8, core.loadCellLine(cellsstring, true, false).getDBID());
						update.setInt(9, core.loadReadType(rtypestring, true, false).getDBID());
						update.setInt(10, readlength);
						update.setInt(11, numreads);
						update.setString(12, collabid);
						update.setString(13, publicsource);
						update.setString(14, publicdbid);
						update.setString(15, fqfile);
						update.setString(16, exptnote);
						update.setInt(17, eID);
						update.execute();
						update.close();
						try {
							expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
						} catch (NotFoundException e2) {
							/* failed again means the insert failed. you lose */
							// cxn.rollback(); //Can't use rollback with
							// autocommit=true?
							throw new DatabaseException(
									"Something went wrong when updating " + alignpieces[0] + "," + alignpieces[1]);
						}
					}
				}

				if (expt == null) {
					// cxn.rollback(); //Can't use rollback with
					// autocommit=true?
					throw new DatabaseException(
							"Couldn't create/update seqexpt " + alignpieces[2] + " for " + alignpieces[0]);
				}
				System.out.println(expt.getDBID());
			} catch (UnknownRoleException e) {
				throw new IllegalArgumentException("Unknown role: seqdata" + e);
			} finally {
				if (insert != null) {
					try {
						insert.close();
					} catch (SQLException ex) {
					}
				}
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
						throw new DatabaseException("Couldn't close connection with role seqdata" + ex);
					}
			}

			loader.close();
		}
	}
}
