package org.seqcode.data.seqdata.tools;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.seqdata.*;
import org.seqcode.utils.*;


/**
 *
 */

public class FindZeroReadAlignments {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        java.sql.Connection cxn = DatabaseConnectionManager.getConnection("seqdata");
        cxn.setAutoCommit(true);
        
        SeqDataLoader loader = new SeqDataLoader();
                
        Collection<SeqExpt> allExpts = loader.loadAllExperiments();
        for(SeqExpt expt : allExpts){
        	if(expt.getNumRead()==0){
        		Collection<SeqAlignment> aligns = loader.loadAlignmentsBySeqExpt(expt);
        		for(SeqAlignment align : aligns){
        			System.out.println(align.getDBID()+"\t"+expt.getName()+"\t"+expt.getReplicate()+"\t"+align.getName()+"\t"+align.getGenome()+"\t"+align.getAlignFile());
        		}
        	}
        }           
        loader.close();
        cxn.close();
    }
}