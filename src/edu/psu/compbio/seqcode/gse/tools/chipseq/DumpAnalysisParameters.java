package edu.psu.compbio.seqcode.gse.tools.chipseq;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class DumpAnalysisParameters {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        ChipSeqAnalysis analysis = Args.parseChipSeqAnalysis(args,"analysis");                                                
        Map<String,String> params = analysis.getParams();
        for (String k : params.keySet()) {
            System.out.println(k + "=" + params.get(k));
        }
    }
}