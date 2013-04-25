package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/** 
 * Creates an experiment (if necessary) and alignment (if necessary)
 * in the database and prints the DBID on stdout.  Use this with
 * the readdb importer to load metadata to the mysql database.
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

	    	java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
	        cxn.setAutoCommit(false);
	        Genome genome = Args.parseGenome(args).cdr();
	        String alignname = Args.parseString(args,"align",null);
	        String alignpieces[] = alignname.split(";");
	        String etypestring = Args.parseString(args,"expttype",null);
	        String labstring = Args.parseString(args,"lab",null);
	        String conditionstring = Args.parseString(args,"exptcondition",null);
	        String targetstring = Args.parseString(args,"expttarget",null);
	        String cellsstring = Args.parseString(args,"cellline",null);
	        String rtypestring = Args.parseString(args,"readtype",null);
	        String atypestring = Args.parseString(args,"aligntype",null);
	        String paramsfname = Args.parseString(args,"paramsfile",null);
	        int readlength = Args.parseInteger(args,"readlength",36);
	        int numreads = Args.parseInteger(args,"numreads",0);
	        String collabid = Args.parseString(args,"collabid",null);
	        String publicsource = Args.parseString(args,"publicsource","UNPUB");
	        String publicdbid = Args.parseString(args,"publicdbid","NA");
	        String fqfile = Args.parseString(args,"fqfile",null);
	        String exptnote = Args.parseString(args,"exptnote",null);
	        String permissions = Args.parseString(args,"permissions",null);
	        int numhits = Args.parseInteger(args,"numhits",0);
	        float totalweight = Args.parseFloat(args,"totalweight",0);
	        int numpairs = Args.parseInteger(args,"numpairs",0);
	        float totalpairweight = Args.parseFloat(args,"totalpairweight",0);
	        String aligndir = Args.parseString(args,"aligndir",null);
	        String alignfile = Args.parseString(args,"alignfile",null);
	        String idxfile = Args.parseString(args,"idxfile",null);
	        String collabalignid = Args.parseString(args,"collabalignid",null);
	        
	        SeqExpt expt = null;
	        SeqAlignment alignment = null;
	        SeqDataLoader loader = new SeqDataLoader();
	        MetadataLoader core = new MetadataLoader();
	        boolean newExpt=false;
	        
	        //SEQEXPERIMENT
	        try {
	            expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
	        } catch (NotFoundException e) {
	        	//NotFound = create new experiment
	        	newExpt=true;
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
	            insert.setString(12, collabid);
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
	        if(!newExpt){
	        	//Experiment exists: Update the old experiment
	        	int eID = expt.getDBID();
	        	System.err.println("Updating experiment "+eID+" " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
	            PreparedStatement insert = SeqExpt.createInsertWithID(cxn);
	            insert.setInt(1, eID);
	            insert.setString(2, alignpieces[0]);
	            insert.setString(3, alignpieces[1]);
	            insert.setInt(4, genome.getSpeciesDBID());
	            insert.setInt(5, core.getExptType(etypestring).getDBID());
	            insert.setInt(6, core.getLab(labstring).getDBID());
	            insert.setInt(7, core.getExptCondition(conditionstring).getDBID());
	            insert.setInt(8, core.getExptTarget(targetstring).getDBID());
	            insert.setInt(9, core.getCellLine(cellsstring).getDBID());
	            insert.setInt(10, core.getReadType(rtypestring).getDBID());
	            insert.setInt(11, readlength);
	            insert.setInt(12, numreads);
	            insert.setString(13, collabid);
	            insert.setString(14, publicsource);
	            insert.setString(15, publicdbid);
	            insert.setString(16, fqfile);
	            insert.setString(17, exptnote);
	            insert.execute();
	            try {
	                expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
	            } catch (NotFoundException e2) {
	                /* failed again means the insert failed.  you lose */
	                cxn.rollback();
	                throw new DatabaseException("Something went wrong when updating " + alignpieces[0] + "," + alignpieces[1]);
	            }
	        }
	        
	        //SEQALIGNMENT
	        alignment = loader.loadAlignment(expt, alignpieces[2], genome);
	        if (alignment == null) {
	            try {
	                PreparedStatement insert = SeqAlignment.createInsertStatement(cxn);
	                System.err.println("Creating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
	                System.err.println("Inserting alignment for experiment " + expt.getDBID());
	                insert.setInt(1, expt.getDBID());
	                insert.setString(2, alignpieces[2]);
	                insert.setInt(3, genome.getDBID());
	                insert.setString(4, permissions);
	                insert.setInt(5, core.getAlignType(atypestring).getDBID());
	                insert.setInt(6, numhits);
	                insert.setFloat(7, totalweight);
	                insert.setInt(8, numpairs);
	                insert.setFloat(9, totalpairweight);
	                insert.setString(10, aligndir);
	                insert.setString(11, alignfile);
	                insert.setString(12, idxfile);
	                insert.setString(13, collabalignid);
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
	        }else{
	        	try {
	        		int aID = alignment.getDBID();
	                PreparedStatement insert = SeqAlignment.createInsertStatementWithID(cxn);
	                System.err.println("Updating alignment "+aID+" " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
	                System.err.println("Updating alignment for experiment " + expt.getDBID());
	                insert.setInt(1, aID);
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
	        
	        
	        core.close();
	        loader.close();
	        if (alignment == null) {
	            cxn.rollback();
	            throw new DatabaseException("Couldn't create/update alignment " + alignpieces[2] + " for " + alignpieces[0]);
	        }
	        System.out.println(alignment.getDBID());
	        cxn.commit();
    	}        
    }
}
