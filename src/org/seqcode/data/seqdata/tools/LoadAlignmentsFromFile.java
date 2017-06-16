package org.seqcode.data.seqdata.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.UnknownRoleException;
import org.seqcode.data.core.MetadataLoader;
import org.seqcode.data.seqdata.SeqAlignment;
import org.seqcode.data.seqdata.SeqDataLoader;
import org.seqcode.data.seqdata.SeqExpt;
import org.seqcode.genome.Genome;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;

/**
 * Load SeqAlignment descriptions from a file. This imports a file like
 * deepseq.list when populating a new installation of the seqdata mysql
 * database.
 * 
 * @author mahony
 *
 *         Usage: LoadAlignmentsFromFile --list "filename"
 * 
 *         The assumed file is in the deepseq.list format, with the following
 *         fields:
 * 
 *         0) ReadDBID 1) ExptType 2) Lab 3) ExptCondition 4) ExptTarget 5)
 *         CellLine 6) Replicate 7) Aligner 8) Genome 9) Permissions 10)
 *         PubSource 11) PublicDBID 12) CollabExptID 13) CollabAlignID 14)
 *         ReadType 15) AlignType 16) ReadLength 17) TotalReads 18) AlignedHits
 *         19) UniquelyAlignedHits 20) DBLoadedHits 21) DBLoadedWeight 22)
 *         DBLoadedType2Hits 23) DBLoadedType2Weight 24) DBLoadedPairs 25)
 *         DBLoadedPairWeight 26) ReadsFile 27) AlignDir 28) AlignFile 29)
 *         IDXFile 30) AlignParamFile 31) ExptNote 32) LoadDate 33) ExptName
 * 
 */
public class LoadAlignmentsFromFile {

	public static void main(String args[]) throws SQLException, IOException, NotFoundException {
		String filename = Args.parseString(args, "list", null);
		if (filename == null) {
			System.err.println("LoadAlignmentsFromFile:\n" + "\t--list <deepseq.list format file>\n");
			System.exit(1);
		}
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));

		SeqDataLoader loader = new SeqDataLoader();
		MetadataLoader core = loader.getMetadataLoader();

		java.sql.Connection cxn = null;
		PreparedStatement insert = null;
		try {
			cxn = DatabaseConnectionManager.getConnection("seqdata");
			cxn.setAutoCommit(true);

			// Iterate through the file
			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if (!fields[0].equals("ReadDBID")) {// skip the first line in
													// the deepseq.list file

					// Variables
					Integer dbid = new Integer(fields[0]);
					String alignpieces[] = fields[33].split(";");
					Genome genome = Genome.findGenome(fields[8]);
					String etypestring = fields[1];
					String labstring = fields[2];
					String conditionstring = fields[3];
					String targetstring = fields[4];
					String cellsstring = fields[5];
					String rtypestring = fields[14];
					String atypestring = fields[15];
					Integer readlength = new Integer(fields[16]);
					String numreadsStr = new String(fields[17]);
					Integer numreads;// Paired reads are counted as one
					if (numreadsStr.contains("+")) {
						String[] tmp = numreadsStr.split("\\+");
						numreads = new Integer(tmp[0]);
					} else
						numreads = new Integer(numreadsStr);
					Integer numhits = new Integer(fields[20]);
					Float totalweight = new Float(fields[21]);
					Integer numtype2hits = new Integer(fields[22]);
					Float totaltype2weight = new Float(fields[23]);
					Integer numpairs = new Integer(fields[24]);
					Float totalpairweight = new Float(fields[25]);
					String permissions = fields[9];
					String collabexptid = fields[12];
					String collabalignid = fields[13];
					String publicsource = fields[10];
					String publicdbid = fields[11];
					String fqfile = fields[26];
					String aligndir = fields[27];
					String alignfile = fields[28];
					String idxfile = fields[29];
					String paramsfname = fields[30];
					String exptnote = fields[31];

					// From here on out, it's similar to CreateAlignment
					SeqExpt expt = null;
					SeqAlignment alignment = null;
					try {
						expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
					} catch (NotFoundException e) {
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
						insert.setString(12, collabexptid);
						insert.setString(13, publicsource);
						insert.setString(14, publicdbid);
						insert.setString(15, fqfile);
						insert.setString(16, exptnote);
						insert.execute();
						try {
							expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
						} catch (NotFoundException e2) {
							/* failed again means the insert failed. you lose */
							// cxn.rollback();
							throw new DatabaseException("Couldn't create " + alignpieces[0] + "," + alignpieces[1]);
						}
					}
					alignment = loader.loadAlignment(expt, alignpieces[2], genome);

					if (alignment == null) {
						try {
							insert = SeqAlignment.createInsertStatementWithID(cxn);
							System.err.println("Creating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";"
									+ alignpieces[2]);
							System.err.println("Inserting alignment for experiment " + expt.getDBID());
							insert.setInt(1, dbid);
							insert.setInt(2, expt.getDBID());
							insert.setString(3, alignpieces[2]);
							insert.setInt(4, genome.getDBID());
							insert.setString(5, permissions);
							insert.setInt(6, core.loadAlignType(atypestring, true, false).getDBID());
							insert.setInt(7, numhits);
							insert.setFloat(8, totalweight);
							insert.setInt(9, numtype2hits);
							insert.setFloat(10, totaltype2weight);
							insert.setInt(11, numpairs);
							insert.setFloat(12, totalpairweight);
							insert.setString(13, aligndir);
							insert.setString(14, alignfile);
							insert.setString(15, idxfile);
							insert.setString(16, collabalignid);
							insert.execute();
							alignment = loader.loadAlignment(expt, alignpieces[2], genome);
							File f = null;
							if (paramsfname != null) {
								f = new File(paramsfname);
							}
							if (f != null && f.exists()) {
								System.err.println("Reading alignment parameters from " + f);
								loader.addAlignmentParameters(alignment, f);

							}
						} catch (IOException e) {
							// cxn.rollback();
							System.err.println("Couldn't add alignment parameters");
							e.printStackTrace();
						}

					}
					if (alignment == null) {
						// cxn.rollback();
						throw new DatabaseException(
								"Couldn't create alignment " + alignpieces[2] + " for " + alignpieces[0]);
					}
					System.out.println(alignment.getDBID());
				}
			}
		} catch (UnknownRoleException e) {
			throw new IllegalArgumentException("Unknown role: seqdata" + e);
		} finally {
			if (insert != null) {
				try {
					insert.close();
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
		reader.close();

	}

}
