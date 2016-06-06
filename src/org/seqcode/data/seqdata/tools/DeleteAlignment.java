package org.seqcode.data.seqdata.tools;

import java.sql.*;
import java.io.IOException;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.core.CellLine;
import org.seqcode.data.core.ExptCondition;
import org.seqcode.data.core.ExptTarget;
import org.seqcode.data.core.Lab;
import org.seqcode.data.seqdata.*;
import org.seqcode.utils.*;


/**
 * Delete the SeqAlignment entry associated with the alignment ID specified on the command line.
 * SeqExpt entries and other tables are edited if the SeqAlignment deletion leaves orphans.
 * 
 * DeleteAlignment --id "999"
 */

public class DeleteAlignment {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
    	java.sql.Connection cxn = DatabaseConnectionManager.getConnection("seqdata");
        cxn.setAutoCommit(true);
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
	        if(loader.loadAlignmentsBySeqExpt(expt).size()==0){
	        	System.err.println("Deleting experiment: "+expt.getName()+"\t"+expt.getDBID());
	        	seqDatamodifier.deleteSeqExpt(expt);
	        }
        }
        loader.close();
        cxn.close();
    }
}