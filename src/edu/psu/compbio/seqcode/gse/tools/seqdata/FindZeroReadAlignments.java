package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 *
 */

public class FindZeroReadAlignments {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
        cxn.setAutoCommit(false);
        
        SeqDataLoader loader = new SeqDataLoader();
                
        Collection<SeqExpt> allExpts = loader.loadAllExperiments();
        for(SeqExpt expt : allExpts){
        	if(expt.getNumRead()==0){
        		Collection<SeqAlignment> aligns = loader.loadAllAlignments(expt);
        		for(SeqAlignment align : aligns){
        			System.out.println(align.getDBID()+"\t"+expt.getName()+"\t"+expt.getReplicate()+"\t"+align.getName()+"\t"+align.getGenome()+"\t"+align.getAlignFile());
        		}
        	}
        }           
        loader.close();
        cxn.close();
    }
}