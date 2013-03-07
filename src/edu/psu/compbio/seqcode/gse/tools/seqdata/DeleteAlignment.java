package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.*;
import java.util.*;
import java.io.IOException;

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
                
        SeqAlignment align = loader.loadAlignment(id);
        SeqExpt expt = align.getExpt();
        
 //TODO       
        //Find and delete Analysis entries (if exist)
        
        //Find and delete the AlignmentParameters (if exist)
        
        //Delete the SeqAlignment
        
        //Delete the SeqExpt if no other SeqAlignments depend
        
        //Delete core.expttype if no other SeqExpts depend
        //Delete core.lab if no other SeqExpts depend
        //Delete core.exptcondition if no other SeqExpts depend
        //Delete core.expttarget if no other SeqExpts depend
        //Delete core.cellline if no other SeqExpts depend
        
             
    }
}