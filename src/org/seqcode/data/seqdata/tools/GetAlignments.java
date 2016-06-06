package org.seqcode.data.seqdata.tools;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import org.seqcode.data.core.*;
import org.seqcode.data.seqdata.*;
import org.seqcode.genome.Genome;
import org.seqcode.utils.*;


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
        
        Genome genome = Args.parseGenome(args).cdr();
        Collection<String> alignnames = Args.parseStrings(args,"align");
        String etypestring = Args.parseString(args,"expttype",null);
        String labstring = Args.parseString(args,"lab",null);
        String targetstring = Args.parseString(args,"expttarget",null);
        String cellsstring = Args.parseString(args,"cellline",null);
        String conditionstring = Args.parseString(args,"exptcondition",null);
        String rtypestring = Args.parseString(args,"readtype",null);

        SeqDataLoader loader = new SeqDataLoader();
        MetadataLoader core = loader.getMetadataLoader();

        Integer target = null, cells = null, condition = null, lab = null, expttype = null, readtype = null;
        if (targetstring != null) {
            target = core.loadExptTarget(targetstring, false, false).getDBID();
        }
        if (cellsstring != null) {
            cells = core.loadCellLine(cellsstring, false, false).getDBID();
        }
        if (conditionstring != null) {
            condition = core.loadExptCondition(conditionstring, false, false).getDBID();
        }
        if (labstring != null) {
            lab = core.loadLab(labstring, false, false).getDBID();
        }
        if (etypestring != null) {
            expttype = core.loadExptType(etypestring, false, false).getDBID();
        }
        if (rtypestring != null) {
            readtype = core.loadReadType(rtypestring, false, false).getDBID();
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
    }
}