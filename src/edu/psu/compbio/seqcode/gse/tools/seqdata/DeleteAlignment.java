package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.general.Lab;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataModifier;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 * Delete the SeqAlignment entry associated with the alignment ID specified on the command line.
 * SeqExpt entries and other tables are edited if the SeqAlignment deletion leaves orphans.
 * 
 * DeleteAlignment --id "999"
 */

public class DeleteAlignment {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
    	java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
        cxn.setAutoCommit(false);
        Integer id = Args.parseInteger(args,"id", -1);
        
        SeqDataLoader loader = new SeqDataLoader();
        SeqDataModifier seqDatamodifier = new SeqDataModifier(loader);
        
        SeqAlignment align = loader.loadAlignment(id);
        
        if(align != null){
	        SeqExpt expt = align.getExpt();
	        Lab lab = expt.getLab();
	        ExptCondition cond = expt.getExptCondition();
	        ExptTarget target = expt.getExptTarget();
	        CellLine cells = expt.getCellLine();
	        
	        //TODO       
	        //Find and delete Analysis entries (if exist)
	        
	        //Find and delete the AlignmentParameters (if exist)
	        System.err.println("Deleting alignment parameters for: "+align.getName());
	        seqDatamodifier.deleteAlignmentParameters(align);
	        
	        //Delete the SeqAlignment
	        System.err.println("Deleting alignment: "+align.getName()+"\t"+align.getDBID());
	        seqDatamodifier.deleteSeqAlignment(align);
	        
	        //Delete the SeqExpt if no other SeqAlignments depend
	        if(loader.loadAllAlignments(expt).size()==0){
	        	System.err.println("Deleting experiment: "+expt.getName()+"\t"+expt.getDBID());
	        	seqDatamodifier.deleteSeqExpt(expt);
	        }
        }
        loader.close();
        cxn.close();
    }
}