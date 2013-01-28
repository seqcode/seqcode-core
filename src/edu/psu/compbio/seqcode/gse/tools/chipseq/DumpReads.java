package edu.psu.compbio.seqcode.gse.tools.chipseq;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 * Generates a FASTA file with reads from experiments that match a specified cells
 * and condition.  
 */

public class DumpReads {


    public static void main(String args[]) throws SQLException, IOException {
        String cellsname = Args.parseString(args,"cells",null);
        String conditionname = Args.parseString(args,"condition",null);
        if (cellsname == null) {
            System.err.println("Must specify --cells");
            System.exit(1);
        }
        if (conditionname == null) {
            System.err.println("Must specify --condition");
            System.exit(1);
        }

        ChipSeqLoader loader = new ChipSeqLoader();
        Collection<ChipSeqExpt> allExpts = loader.loadAllExperiments();

        Connection cxn = DatabaseFactory.getConnection("chipseq");
        PreparedStatement stmt = cxn.prepareStatement("select name, sequence from chipseqreads where expt = ?");

        for (ChipSeqExpt expt : allExpts) {
            System.err.println("Considering " + expt);
            if (expt.getCells().getName().equals(cellsname) &&
                expt.getCondition().getName().equals(conditionname)) {
                stmt.setInt(1, expt.getDBID());
                ResultSet rs = stmt.executeQuery();
                System.err.println("dumping");
                while (rs.next()) {
                    System.out.println(">" + rs.getString(1) + "\n" + rs.getString(2));
                }
                rs.close();
            }
        }
        stmt.close();        
        loader.close();
    }
}
