package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 * Returns the SeqExpt and SeqAlignment data associated with the alignment IDs specified on the command line.
 * 
 * DescribeAlignments --id "999"
 */

public class DescribeAlignments {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
        cxn.setAutoCommit(false);
        Collection<Integer> ids = Args.parseIntegers(args,"id");
        
        SeqDataLoader loader = new SeqDataLoader();
                
        for (Integer id : ids) {
            SeqAlignment align = loader.loadAlignment(id);
            SeqExpt expt = align.getExpt();
            
            System.out.println("----------\nID="+id);
            
            System.out.println("SeqExpt info:");
            System.out.println("\tID: "+expt.getDBID());
            System.out.println("\tName: "+expt.getName());
            System.out.println("\tRep: "+expt.getReplicate());
            System.out.println("\tSpecies: "+expt.getOrganism().getName());
            System.out.println("\tLab: "+expt.getLab().getName());
            System.out.println("\tExptType: "+expt.getExptType().getName());
            System.out.println("\tExptCondition: "+expt.getExptCondition().getName());
            System.out.println("\tExptTarget: "+expt.getExptTarget().getName());
            System.out.println("\tExptCellLine: "+expt.getCellLine().getName());
            System.out.println("\tReadType: "+expt.getReadType().getName());
            System.out.println("\tReadLen: "+expt.getReadLength());
            System.out.println("\tNumRead: "+expt.getNumRead());
            System.out.println("\tCollabID: "+expt.getCollabID());
            System.out.println("\tPublicSource: "+expt.getPublicSource());
            System.out.println("\tPublicDBID: "+expt.getPublicDBID());
            System.out.println("\tFQFile: "+expt.getFQFile());
            System.out.println("\tExptNote: "+expt.getExptNote());
            
            System.out.println("SeqAlignment info:");
            System.out.println("\tID: "+align.getDBID());
            System.out.println("\tName: "+align.getName());
            System.out.println("\tGenome: "+align.getGenome().getVersion());
            System.out.println("\tAlignType: "+align.getAlignType().getName());
            System.out.println("\tPermissions: "+align.getPermissions());
            System.out.println("\tNumHits: "+align.getNumHits());
            System.out.println("\tTotalWeight: "+align.getTotalWeight());
            System.out.println("\tAlignDir: "+align.getAlignDir());
            System.out.println("\tAlignFile: "+align.getAlignFile());
            System.out.println("\tIDXFile: "+align.getIDXFile());
            System.out.println("\tCollabAlignID: "+align.getCollabAlignID());
            
        }
        cxn.close();
    }
}