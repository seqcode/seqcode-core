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
public class UpdateAlignmentsFromFile {

    public static void main(String args[]) throws SQLException, IOException, NotFoundException {
    	String filename = Args.parseString(args,"list",null);
    	if(filename==null){
    		System.err.println("UpdateAlignmentsFromFile:\n" +
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
			if(!fields[0].equals("ReadDBID") && !fields[0].startsWith("#")){//skip the first line in the deepseq.list file
				
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
				
				//Have to load alignment and experiment by DBID since the naming may change
				SeqExpt expt = null;
		        SeqAlignment alignment = loader.loadAlignment(dbid);
		        if (alignment == null) {
		        	cxn.rollback();
		            throw new DatabaseException("Can't find alignment "+dbid+" for " + alignpieces[2] + " for " + alignpieces[0]);
		        }else{
		        	//Update experiment
		        	int exptID = alignment.getExpt().getDBID();
			        boolean exptExists=false;
			        try {
			            expt = loader.loadExperiment(exptID);
			            exptExists=true;
			        } catch (NotFoundException e) {
			        	cxn.rollback();
			            System.err.println("No experiment found for " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
			            System.exit(1);
			        }
			        if(exptExists){
			        	System.err.println("Updating experiment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
			            PreparedStatement update = SeqExpt.createUpdateWithID(cxn);
			            update.setString(1, alignpieces[0]);
			            update.setString(2, alignpieces[1]);
			            update.setInt(3, genome.getSpeciesDBID());
			            update.setInt(4, core.getExptType(etypestring).getDBID());
			            update.setInt(5, core.getLab(labstring).getDBID());
			            update.setInt(6, core.getExptCondition(conditionstring).getDBID());
			            update.setInt(7, core.getExptTarget(targetstring).getDBID());
			            update.setInt(8, core.getCellLine(cellsstring).getDBID());
			            update.setInt(9, core.getReadType(rtypestring).getDBID());
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
			                /* failed again means the insert failed.  you lose */
			                cxn.rollback();
			                throw new DatabaseException("Couldn't update experiment for " + alignpieces[0] + "," + alignpieces[1]);
			            }
			        }
			       
			        
			        //Alignment already loaded above
		            try {
		                PreparedStatement update = SeqAlignment.createUpdateStatementWithID(cxn);
		                System.err.println("Updating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
		                update.setInt(1, expt.getDBID());
		                update.setString(2, alignpieces[2]);
		                update.setInt(3, genome.getDBID());
		                update.setString(4, permissions);
		                update.setInt(5, core.getAlignType(atypestring).getDBID());
		                update.setInt(6, numhits);
		                update.setFloat(7, totalweight);
		                update.setInt(8, numpairs);
		                update.setFloat(9, totalpairweight);
		                update.setString(10, aligndir);
		                update.setString(11, alignfile);
		                update.setString(12, idxfile);
		                update.setString(13, collabalignid);
		                update.setInt(14, dbid);
		                update.execute();
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

