package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class DumpAnalysisParameters {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        SeqAnalysis analysis = Args.parseSeqAnalysis(args,"analysis");                                                
        Map<String,String> params = analysis.getParams();
        for (String k : params.keySet()) {
            System.out.println(k + "=" + params.get(k));
        }
    }
}