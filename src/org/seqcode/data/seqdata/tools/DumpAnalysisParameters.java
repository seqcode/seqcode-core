package org.seqcode.data.seqdata.tools;

import java.util.*;

import org.seqcode.data.seqdata.*;
import org.seqcode.genome.Genome;
import org.seqcode.gseutils.Args;


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