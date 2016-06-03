package org.seqcode.gse.tools.seqdata;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import org.seqcode.genome.Genome;
import org.seqcode.gse.datasets.core.MetadataLoader;
import org.seqcode.gse.datasets.seqdata.SeqAlignment;
import org.seqcode.gse.datasets.seqdata.SeqDataLoader;
import org.seqcode.gse.datasets.seqdata.SeqExpt;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.database.DatabaseConnectionManager;
import org.seqcode.gse.utils.database.DatabaseException;
import org.seqcode.gse.utils.database.UnknownRoleException;


/**
 * Update existing SeqExpt & SeqAlignment descriptions from a file. 
 * This imports a file like deepseq.list when populating a new installation 
 * of the seqdata mysql database. 
 *  
 * @author mahony
 *
 * Usage: UpdateAlignmentsFromFile --list "filename"
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
 *22) DBLoadedType2Hits
 *23) DBLoadedType2Weight
 *24) DBLoadedPairs
 *25) DBLoadedPairWeight
 *26) ReadsFile
 *27) AlignDir
 *28) AlignFile
 *29) IDXFile
 *30) AlignParamFile
 *31) ExptNote
 *32) LoadDate
 *33) ExptName
 * 
 */
public class UpdateAlignmentsFromFile {

    public static void main(String args[]) throws SQLException, IOException, NotFoundException {
    	String filename = Args.parseString(args,"list",null);
    	if(filename==null){
    		System.err.println("UpdateAlignmentsFromFile:\n" +
    				"\t--list <deepseq.list format file>\n");
    		System.exit(1);
    	}
    	BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
        
        SeqDataLoader loader = new SeqDataLoader();
        MetadataLoader core = loader.getMetadataLoader();

        java.sql.Connection cxn = null;
        PreparedStatement update = null;
        try{
        	cxn = DatabaseConnectionManager.getConnection("seqdata");
            cxn.setAutoCommit(true);
            
	        //Iterate through the file
	        String line = null;
			while ((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if(!fields[0].equals("ReadDBID") && !fields[0].startsWith("#")){//skip the first line in the deepseq.list file
					
					//Variables
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
					Integer numreads;//Paired reads are counted as one
					if(numreadsStr.contains("+")){
						String[] tmp = numreadsStr.split("\\+");
						numreads = new Integer(tmp[0]);
					}else
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
					
					//Have to load alignment and experiment by DBID since the naming may change
					SeqExpt expt = null;
			        SeqAlignment alignment = loader.loadAlignment(dbid);
			        if (alignment == null) {
			        	//cxn.rollback();
			            throw new DatabaseException("Can't find alignment "+dbid+" for " + alignpieces[2] + " for " + alignpieces[0]);
			        }else{
			        	//Update experiment
			        	int exptID = alignment.getExpt().getDBID();
				        boolean exptExists=false;
				        try {
				            expt = loader.loadExperiment(exptID);
				            exptExists=true;
				        } catch (NotFoundException e) {
				        	//cxn.rollback();
				            System.err.println("No experiment found for " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
				            System.exit(1);
				        }
				        if(exptExists){
				        	System.err.println("Updating experiment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
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
				            update.setString(12, collabexptid);
				            update.setString(13, publicsource);
				            update.setString(14, publicdbid);
				            update.setString(15, fqfile);
				            update.setString(16, exptnote);
				            update.setInt(17, expt.getDBID());
				            update.execute();
				            try {
				                expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
				            } catch (NotFoundException e2) {
				                /* failed again means the update failed.  you lose */
				                //cxn.rollback();
				                throw new DatabaseException("Couldn't update experiment for " + alignpieces[0] + "," + alignpieces[1]);
				            }
				        }
				       
				        
				        //Alignment already loaded above
			            try {
			                update = SeqAlignment.createUpdateStatementWithID(cxn);
			                System.err.println("Updating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
			                update.setInt(1, expt.getDBID());
			                update.setString(2, alignpieces[2]);
			                update.setInt(3, genome.getDBID());
			                update.setString(4, permissions);
			                update.setInt(5, core.loadAlignType(atypestring, true, false).getDBID());
			                update.setInt(6, numhits);
			                update.setFloat(7, totalweight);
			                update.setInt(8, numtype2hits);
			                update.setFloat(9, totaltype2weight);
			                update.setInt(10, numpairs);
			                update.setFloat(11, totalpairweight);
			                update.setString(12, aligndir);
			                update.setString(13, alignfile);
			                update.setString(14, idxfile);
			                update.setString(15, collabalignid);
			                update.setInt(16, dbid);
			                update.execute();
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
			                //cxn.rollback();
			                System.err.println("Couldn't add alignment parameters");
			                e.printStackTrace();
			            }
	
			        }
			        System.out.println(alignment.getDBID());
				}
			}
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: seqdata" + e);
        } finally {
        	if (update != null) { try {update.close(); } catch (SQLException ex) {  }}
	        if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role seqdata"+ ex); }
        }
		loader.close();
		reader.close();
    }
    
}

