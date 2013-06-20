package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;

/**
 * Load SeqAlignment descriptions from a file. 
 * This imports a file like deepseq.list when populating a new installation 
 * of the seqdata mysql database. 
 *  
 * @author mahony
 *
 * Usage: LoadAlignmentsFromFile --list "filename"
 * 
 * The assumed file is in the deepseq.list format, with the following fields:
 * 
 *0) ReadDBID
 *1) ExptType
 *2) Lab
 *3) ExptCondition
 *4) ExptTarget
 *5) CellLine
 *6) Replicate
 *7) Aligner
 *8) Genome
 *9) Permissions
 *10) PubSource
 *11) PublicDBID
 *12) CollabExptID
 *13) CollabAlignID
 *14) ReadType
 *15) AlignType
 *16) ReadLength
 *17) TotalReads
 *18) AlignedHits
 *19) UniquelyAlignedHits
 *20) DBLoadedHits
 *21) DBLoadedWeight
 *22) DBLoadedPairs
 *23) DBLoadedPairWeight
 *24) ReadsFile
 *25) AlignDir
 *26) AlignFile
 *27) IDXFile
 *28) AlignParamFile
 *29) ExptNote
 *30) LoadDate
 *31) ExptName
 * 
 */
public class LoadAlignmentsFromFile {

    public static void main(String args[]) throws SQLException, IOException, NotFoundException {
    	String filename = Args.parseString(args,"list",null);
    	if(filename==null){
    		System.err.println("LoadAlignmentsFromFile:\n" +
    				"\t--list <deepseq.list format file>\n");
    		System.exit(1);
    	}
        java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
        cxn.setAutoCommit(false);
        
        SeqDataLoader loader = new SeqDataLoader();
        MetadataLoader core = new MetadataLoader();

        //Iterate through the file
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
        String line = null;
		while ((line = reader.readLine()) != null) {
			String[] fields = line.split("\t");
			if(!fields[0].equals("ReadDBID")){//skip the first line in the deepseq.list file
				
				//Variables
				Integer dbid = new Integer(fields[0]);
				String alignpieces[] = fields[31].split(";");
				Genome genome = Organism.findGenome(fields[8]);
				String etypestring = fields[1];
				String labstring = fields[2];
				String conditionstring = fields[3];
				String targetstring = fields[4];
				String cellsstring = fields[5];
				String rtypestring = fields[14];
				String atypestring = fields[15];
				Integer readlength = new Integer(fields[16]);
				String numreadsStr = new String(fields[17]);
				Integer numreads;//Paired reads are counted as one
				if(numreadsStr.contains("+")){
					String[] tmp = numreadsStr.split("\\+");
					numreads = new Integer(tmp[0]);
				}else
					numreads = new Integer(numreadsStr);
				Integer numhits = new Integer(fields[20]);
				Float totalweight = new Float(fields[21]);
				Integer numpairs = new Integer(fields[22]);
				Float totalpairweight = new Float(fields[23]);
				String permissions = fields[9];
				String collabexptid = fields[12];
				String collabalignid = fields[13];
				String publicsource = fields[10];
				String publicdbid = fields[11];
				String fqfile = fields[24];
				String aligndir = fields[25];
				String alignfile = fields[26];
				String idxfile = fields[27];
				String paramsfname = fields[28];
				String exptnote = fields[29];
				
				//From here on out, it's similar to CreateAlignment
				SeqExpt expt = null;
		        SeqAlignment alignment = null;
		        try {
		            expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
		        } catch (NotFoundException e) {
		            System.err.println("Creating experiment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
		            PreparedStatement insert = SeqExpt.createInsert(cxn);
		            insert.setString(1, alignpieces[0]);
		            insert.setString(2, alignpieces[1]);
		            insert.setInt(3, genome.getSpeciesDBID());
		            insert.setInt(4, core.getExptType(etypestring).getDBID());
		            insert.setInt(5, core.getLab(labstring).getDBID());
		            insert.setInt(6, core.getExptCondition(conditionstring).getDBID());
		            insert.setInt(7, core.getExptTarget(targetstring).getDBID());
		            insert.setInt(8, core.getCellLine(cellsstring).getDBID());
		            insert.setInt(9, core.getReadType(rtypestring).getDBID());
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
		                /* failed again means the insert failed.  you lose */
		                cxn.rollback();
		                throw new DatabaseException("Couldn't create " + alignpieces[0] + "," + alignpieces[1]);
		            }
		        }
		        alignment = loader.loadAlignment(expt, alignpieces[2], genome);

		        if (alignment == null) {
		            try {
		                PreparedStatement insert = SeqAlignment.createInsertStatementWithID(cxn);
		                System.err.println("Creating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
		                System.err.println("Inserting alignment for experiment " + expt.getDBID());
		                insert.setInt(1, dbid);
		                insert.setInt(2, expt.getDBID());
		                insert.setString(3, alignpieces[2]);
		                insert.setInt(4, genome.getDBID());
		                insert.setString(5, permissions);
		                insert.setInt(6, core.getAlignType(atypestring).getDBID());
		                insert.setInt(7, numhits);
		                insert.setFloat(8, totalweight);
		                insert.setInt(9, numpairs);
		                insert.setFloat(10, totalpairweight);
		                insert.setString(11, aligndir);
		                insert.setString(12, alignfile);
		                insert.setString(13, idxfile);
		                insert.setString(14, collabalignid);
		                insert.execute();
		                alignment = loader.loadAlignment(expt, alignpieces[2], genome);
		                cxn.commit();
		                File f = null;
		                if (paramsfname != null) {
		                    f = new File(paramsfname);
		                }
		                if (f != null && f.exists()) {
		                    System.err.println("Reading alignment parameters from " + f);
		                    loader.addAlignmentParameters(alignment, f);

		                }
					} catch (IOException e) {
		                cxn.rollback();
		                System.err.println("Couldn't add alignment parameters");
		                e.printStackTrace();
		            }

		        }
		        if (alignment == null) {
		            cxn.rollback();
		            throw new DatabaseException("Couldn't create alignment " + alignpieces[2] + " for " + alignpieces[0]);
		        }
		        System.out.println(alignment.getDBID());
		        cxn.commit();
			}
		}
		loader.close();
        core.close();
		reader.close();
		cxn.close();
    }
    
}

