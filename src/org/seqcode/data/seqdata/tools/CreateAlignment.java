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


/** 
 * Creates an alignment (if necessary)
 * in the database and prints the DBID on stdout.  Use this with
 * the readdb importer to load metadata to the mysql database.
 * I split this out from a combined experiment/alignment creator because of what seems like
 * simultaneous write conflicts during experiment creation.
 * 
 * Usage:
 * CreateAlignment --species "$SC;SGDv1" --align "name;replicate;alignment version" --lab "Pugh" --expttype "CHIPSEQ" --expttarget "Gcn4" --cellline "FY4" --exptcondition "YPD" --readtype "SINGLE" --aligntype "SINGLE" --paramsfile params.txt --readlength 36
 */
public class CreateAlignment {

    public static void main(String args[]) throws SQLException, IOException, NotFoundException {

    	if(args.length==0){
    		System.out.println("CreateAlignment:\n" +
    				"\t--species <species;genome>\n" +
    				"\t--align <name;replicate;version>\n" +
    				"\t--expttype <CHIPSEQ/CHIPEXO/RNASEQ/etc>\n" +
    				"\t--lab <name>\n" +
    				"\t--exptcondition <condition>\n" +
    				"\t--expttarget <target>\n" +
    				"\t--cellline <cell line>\n" +
    				"\t--readtype <SINGLE/PAIRED>\n" +
    				"\t--aligntype <SINGLE/PAIRED>\n" +
    				"\t--paramsfile <filename>\n" +
    				"\t--readlength <int>\n" +
    				"\t--numreads <int>\n" +
    				"\t--numhits <int>\n" +
    				"\t--totalweight <float>\n" +
    				"\t--numtype2hits <int>\n" +
    				"\t--totaltype2weight <float>\n" +
    				"\t--numpairs <int>\n" +
    				"\t--totalpairweight <float>\n" +
    				"\t--collabid <expt ID>\n" +
    				"\t--collabalignid <align ID>\n" +
    				"\t--publicsource <PMID/UNPUB>\n" +
    				"\t--publicdbid <GEO ID>\n" +
    				"\t--fqfile <FQ filename>\n" +
    				"\t--exptnote <notes about expt>\n" +
    				"\t--permissions <mahony;mahonylab;etc>\n" +
    				"\t--aligndir <directory name>\n" +
    				"\t--alignfile <file name>\n" +
    				"\t--idxfile <file name>\n");
    	}else{

    		java.sql.Connection cxn = null;
    		PreparedStatement insert = null;
    		PreparedStatement update = null;
    		
    		SeqExpt expt = null;
	        SeqAlignment alignment = null;
	        SeqDataLoader loader = new SeqDataLoader();
	        MetadataLoader core = loader.getMetadataLoader();
	        
    		try{
		    	cxn = DatabaseConnectionManager.getConnection("seqdata");
		        cxn.setAutoCommit(true);
		        Genome genome = Args.parseGenome(args).cdr();
		        String alignname = Args.parseString(args,"align",null);
		        String alignpieces[] = alignname.split(";");
		        String labstring = Args.parseString(args,"lab",null);
		        String conditionstring = Args.parseString(args,"exptcondition",null);
		        String targetstring = Args.parseString(args,"expttarget",null);
		        String cellsstring = Args.parseString(args,"cellline",null);
		        String rtypestring = Args.parseString(args,"readtype",null);
		        String atypestring = Args.parseString(args,"aligntype",null);
		        String paramsfname = Args.parseString(args,"paramsfile",null);
		        String permissions = Args.parseString(args,"permissions",null);
		        int numhits = Args.parseInteger(args,"numhits",0);
		        float totalweight = Args.parseFloat(args,"totalweight",0);
		        int numtype2hits = Args.parseInteger(args,"numtype2hits",0);
		        float totaltype2weight = Args.parseFloat(args,"totaltype2weight",0);
		        int numpairs = Args.parseInteger(args,"numpairs",0);
		        float totalpairweight = Args.parseFloat(args,"totalpairweight",0);
		        String aligndir = Args.parseString(args,"aligndir",null);
		        String alignfile = Args.parseString(args,"alignfile",null);
		        String idxfile = Args.parseString(args,"idxfile",null);
		        String collabalignid = Args.parseString(args,"collabalignid",null);
		        
		        int aligntypeID = core.loadAlignType(atypestring, true, false).getDBID();
		        
		        //LOAD THE SEQEXPERIMENT
		        try {
		            expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
		        } catch (NotFoundException e) {
		        	System.err.println("No such experiment: use CreateExpt first");
		        	System.exit(1);
		        }
		        
		        
		        //CREATE THE SEQALIGNMENT
		        alignment = loader.loadAlignment(expt, alignpieces[2], genome);
		        if (alignment == null) {
		            try {
		                insert = SeqAlignment.createInsertStatement(cxn);
		                System.err.println("Creating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
		                System.err.println("Inserting alignment for experiment " + expt.getDBID());
		                insert.setInt(1, expt.getDBID());
		                insert.setString(2, alignpieces[2]);
		                insert.setInt(3, genome.getDBID());
		                insert.setString(4, permissions);
		                insert.setInt(5, core.loadAlignType(atypestring, true, false).getDBID());
		                insert.setInt(6, numhits);
		                insert.setFloat(7, totalweight);
		                insert.setInt(8, numtype2hits);
		                insert.setFloat(9, totaltype2weight);
		                insert.setInt(10, numpairs);
		                insert.setFloat(11, totalpairweight);
		                insert.setString(12, aligndir);
		                insert.setString(13, alignfile);
		                insert.setString(14, idxfile);
		                insert.setString(15, collabalignid);
		                insert.execute();
		                insert.close();
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
		        }else{
		        	//Check if the alignment information is the same before bothering to update the database. 
		        	if(!(alignment.getName().equals(alignpieces[2]) && alignment.getPermissions().equals(permissions) &&
		        		alignment.getAlignType().getDBID()==aligntypeID && alignment.getNumHits()==numhits &&
		        		alignment.getTotalWeight()==totalweight && alignment.getNumPairs()==numpairs && alignment.getTotalPairWeight()==totalpairweight &&
		        		alignment.getAlignDir().equals(aligndir) && alignment.getAlignFile().equals(alignfile) && 
		        		alignment.getIDXFile().equals(idxfile) && alignment.getCollabAlignID().equals(collabalignid))){
			        	try {
			        		int aID = alignment.getDBID();
			                update = SeqAlignment.createUpdateStatementWithID(cxn);
			                System.err.println("Updating alignment "+aID+" " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
			                System.err.println("Updating alignment for experiment " + expt.getDBID());
			                update.setInt(1, expt.getDBID());
			                update.setString(2, alignpieces[2]);
			                update.setInt(3, genome.getDBID());
			                update.setString(4, permissions);
			                update.setInt(5, aligntypeID);
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
			                update.setInt(16, aID);
			                update.execute();
			                update.close();
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
		        }
		        if (alignment == null) {
		            //cxn.rollback();
		            throw new DatabaseException("Couldn't create/update alignment " + alignpieces[2] + " for " + alignpieces[0]);
		        }
		        System.out.println(alignment.getDBID());
    		} catch (UnknownRoleException e) {
                throw new IllegalArgumentException("Unknown role: seqdata" + e);
            } finally {
            	if (insert != null) { try {insert.close(); } catch (SQLException ex) {  }}
    	        if (update != null) { try {update.close();} catch (SQLException ex) { } }
            	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role seqdata"+ ex); }
            }

    		loader.close();
    	}        
    }
}
