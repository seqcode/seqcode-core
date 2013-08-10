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
 * Creates an experiment (if necessary) in the database and prints the DBID on stdout.  
 * Use this with the readdb importer to load metadata to the mysql database.
 * I split this out from a combined experiment/alignment creator because of what seems like
 * simultaneous write conflicts during experiment creation.
 *
 * Usage:
 * CreateExpt --species "$SC;SGDv1" --align "name;replicate;alignment version" --lab "Pugh" --expttype "CHIPSEQ" --expttarget "Gcn4" --cellline "FY4" --exptcondition "YPD" --readtype "SINGLE" --aligntype "SINGLE" --paramsfile params.txt --readlength 36
 */
public class CreateExpt {

    public static void main(String args[]) throws SQLException, IOException, NotFoundException {

    	if(args.length==0){
    		System.out.println("CreateExpt:\n" +
    				"\t--species <species;genome>\n" +
    				"\t--align <name;replicate;version>\n" +
    				"\t--expttype <CHIPSEQ/CHIPEXO/RNASEQ/etc>\n" +
    				"\t--lab <name>\n" +
    				"\t--exptcondition <condition>\n" +
    				"\t--expttarget <target>\n" +
    				"\t--cellline <cell line>\n" +
    				"\t--readtype <SINGLE/PAIRED>\n" +
    				"\t--readlength <int>\n" +
    				"\t--numreads <int>\n" +
    				"\t--collabid <expt ID>\n" +
    				"\t--publicsource <PMID/UNPUB>\n" +
    				"\t--publicdbid <GEO ID>\n" +
    				"\t--fqfile <FQ filename>\n" +
    				"\t--exptnote <notes about expt>\n");
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
	        int readlength = Args.parseInteger(args,"readlength",36);
	        int numreads = Args.parseInteger(args,"numreads",0);
	        String collabid = Args.parseString(args,"collabid",null);
	        String publicsource = Args.parseString(args,"publicsource","UNPUB");
	        String publicdbid = Args.parseString(args,"publicdbid","NA");
	        String fqfile = Args.parseString(args,"fqfile",null);
	        String exptnote = Args.parseString(args,"exptnote",null);
	        
	        SeqExpt expt = null;
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
	            insert.close();
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
	        	//Check if the experiment information is the same before bothering to update the database.
	        	if(!(expt.getName().equals(alignpieces[0]) && expt.getReplicate().equals(alignpieces[1]) &&
	        		expt.getExptType().getDBID()==core.getExptType(etypestring).getDBID() && 
	        		expt.getLab().getDBID()==core.getLab(labstring).getDBID() &&
	        		expt.getExptCondition().getDBID()==core.getExptCondition(conditionstring).getDBID() &&
	        		expt.getExptTarget().getDBID()==core.getExptTarget(targetstring).getDBID() && 
	        		expt.getCellLine().getDBID()==core.getCellLine(cellsstring).getDBID() && 
	        		expt.getReadType().getDBID()==core.getReadType(rtypestring).getDBID() &&
	        		expt.getReadLength()==readlength && expt.getNumRead()==numreads && 
	        		expt.getCollabID().equals(collabid) && expt.getPublicSource().equals(publicsource) &&
	        		expt.getPublicDBID().equals(publicdbid) && expt.getFQFile().equals(fqfile) &&
	        		expt.getExptNote().equals(exptnote))){
	        		
		        	int eID = expt.getDBID();
		        	System.err.println("Updating experiment "+eID+" " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
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
		                /* failed again means the insert failed.  you lose */
		                cxn.rollback();
		                throw new DatabaseException("Something went wrong when updating " + alignpieces[0] + "," + alignpieces[1]);
		            }
	        	}
	        }
	        
	        
	        core.close();
	        loader.close();
	        if (expt == null) {
	            cxn.rollback();
	            throw new DatabaseException("Couldn't create/update seqexpt " + alignpieces[2] + " for " + alignpieces[0]);
	        }
	        System.out.println(expt.getDBID());
	        genome.close();
	        cxn.commit();
	        cxn.close();
    	}        
    }
}
