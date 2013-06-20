package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.general.Lab;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataDeleter;
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
        MetadataDeleter metaDeleter = new MetadataDeleter();
        
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
	        loader.deleteAlignmentParameters(align);
	        
	        //Delete the SeqAlignment
	        System.err.println("Deleting alignment: "+align.getName()+"\t"+align.getDBID());
	        PreparedStatement deleteAlign = SeqAlignment.createDeleteByIDStatement(cxn);
	        deleteAlign.setInt(1, align.getDBID());
	        deleteAlign.execute();
	        deleteAlign.close();
	        cxn.commit();
	        
	        //Delete the SeqExpt if no other SeqAlignments depend
	        if(loader.loadAllAlignments(expt).size()==0){
	        	System.err.println("Deleting experiment: "+expt.getName()+"\t"+expt.getDBID());
	        	PreparedStatement deleteExpt = SeqExpt.createDeleteByDBID(cxn);
	        	deleteExpt.setInt(1, expt.getDBID());
	        	deleteExpt.execute();
	        	deleteExpt.close();
	        	cxn.commit();
	        	
	        	//Delete core.lab if no other SeqExpts depend
	        	if(loader.loadExperiments(lab).size()==0){
	        		System.err.println("Deleting lab: "+lab.getName()+"\t"+lab.getDBID());
	        		metaDeleter.deleteLab(lab.getDBID());
	        	}
	        	
		        //Delete core.exptcondition if no other SeqExpts depend
	        	if(loader.loadExperiments(cond).size()==0){
	        		System.err.println("Deleting condition: "+cond.getName()+"\t"+cond.getDBID());
	        		metaDeleter.deleteCond(cond.getDBID());
	        	}
	        	
		        //Delete core.expttarget if no other SeqExpts depend
	        	if(loader.loadExperiments(target).size()==0){
	        		System.err.println("Deleting target: "+target.getName()+"\t"+target.getDBID());
	        		metaDeleter.deleteTarget(target.getDBID());
	        	}
	        	
		        //Delete core.cellline if no other SeqExpts depend
	        	if(loader.loadExperiments(cells).size()==0){
	        		System.err.println("Deleting cell-line: "+cells.getName()+"\t"+cells.getDBID());
	        		metaDeleter.deleteCell(cells.getDBID());
	        	}
	        }
        }
        loader.close();
        metaDeleter.close();
        cxn.close();
    }
}