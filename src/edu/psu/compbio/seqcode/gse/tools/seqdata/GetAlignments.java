package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 * Returns the ids of the alignments specified on the command line.
 * All fields are optional.  expttype, lab, expttarget, cellline, exptcondition, readtype are ANDed together.  
 * align parameters are ORed together and then applied to the expttype/lab/expttarget/cellline/exptcondition/readtype.  The
 * fields include --align name;replicate;version
 *                --expttype
 *                --lab
 *                --expttarget
 *                --cellline
 *                --exptcondition
 *                --readtype
 *
 * GetAlignments --species "$SC;SGDv1" --align "name;replicate;alignment version" --expttarget "Gcn4" --exptcondition "YPD"
 */

public class GetAlignments {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
        cxn.setAutoCommit(false);
        Genome genome = Args.parseGenome(args).cdr();
        Collection<String> alignnames = Args.parseStrings(args,"align");
        String etypestring = Args.parseString(args,"expttype",null);
        String labstring = Args.parseString(args,"lab",null);
        String targetstring = Args.parseString(args,"expttarget",null);
        String cellsstring = Args.parseString(args,"cellline",null);
        String conditionstring = Args.parseString(args,"exptcondition",null);
        String rtypestring = Args.parseString(args,"readtype",null);

        SeqDataLoader loader = new SeqDataLoader();
        MetadataLoader core = new MetadataLoader();

        Integer target = null, cells = null, condition = null, lab = null, expttype = null, readtype = null;
        if (targetstring != null) {
            target = core.getExptTarget(targetstring).getDBID();
        }
        if (cellsstring != null) {
            cells = core.getCellLine(cellsstring).getDBID();
        }
        if (conditionstring != null) {
            condition = core.getExptCondition(conditionstring).getDBID();
        }
        if (labstring != null) {
            lab = core.getLab(labstring).getDBID();
        }
        if (etypestring != null) {
            expttype = core.getExptType(etypestring).getDBID();
        }
        if (rtypestring != null) {
            readtype = core.getReadType(rtypestring).getDBID();
        }

        for (String an : alignnames) {
            String pieces[] = an.split(";");
            for (SeqAlignment a : loader.loadAlignments(pieces[0],
                                                            (pieces.length > 1 && pieces[1] != "") ? pieces[1] : null,
                                                            (pieces.length > 2 && pieces[2] != "") ? pieces[2] : null,
                                                            expttype, lab, condition,
                                                            target,cells,readtype,genome)) {
                System.out.println(a.getDBID());
            }
        }
        
        loader.close();
        core.close();
        genome.close();
        cxn.close();
    }
}