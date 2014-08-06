package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.util.*;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;

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