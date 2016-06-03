package org.seqcode.gse.tools.seqdata;

import java.util.*;

import org.seqcode.genome.Genome;
import org.seqcode.gse.datasets.seqdata.*;
import org.seqcode.gse.tools.utils.Args;


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