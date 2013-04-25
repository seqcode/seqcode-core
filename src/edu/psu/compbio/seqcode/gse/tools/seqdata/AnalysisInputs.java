package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.util.Collection;
import java.util.Collections;
import java.util.ArrayList;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;


public class AnalysisInputs {
    
    public static void main(String[] args) throws Exception {
        Collection<SeqAnalysis> analyses = Args.parseSeqAnalyses(args,"chipseqanalysis");
        for (SeqAnalysis analysis : analyses) {
            for (SeqAlignment align : analysis.getForeground()) {
                SeqExpt expt = align.getExpt();
                System.out.println(analysis.toString() + "\t" + "FOREGROUND" + "\t" + 
                                   expt.getName() + ";" + expt.getReplicate() + ";" + align.getName());
            }
            for (SeqAlignment align : analysis.getBackground()) {
                SeqExpt expt = align.getExpt();
                System.out.println(analysis.toString() + "\t" + "BACKGROUND" + "\t" + 
                                   expt.getName() + ";" + expt.getReplicate() + ";" + align.getName());
            }

        }

    }


}