package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.util.*;
import java.sql.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;



/**
 * java edu.psu.compbio.seqcode.gse.tools.chipseq.GeneBasedBindingReport --species "$MM;mm9" \
 * --analysisname "PPG ES iCdx2 p2A 7-28-10 lane 5 (36bp)"              \
 * --analysisversion "vs PPG Day4 null-antiV5 iTF_iOlig2 1 (default params) run 2 round 3" \
 * --genes refGene [--proxup 5000] [--proxdown200] [--up 10000] [--intronlen 10000] [--thresh .001]
 *
 * Output columns are
 * 0) gene name
 * 1) positions of distal binding events
 * 2) positions of proximal binding events
 * 3) positions of intronic or exonic binding events
 */


public class GeneBasedBindingReport extends GeneBasedReport {

    private SeqAnalysis analysis;
    private double pvalthresh;

    public static void main(String args[]) throws Exception {
        GeneBasedBindingReport report = new GeneBasedBindingReport();
        report.parseArgs(args);
        report.getRegions(args);
        report.report();
    }

    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        super.parseArgs(args);
        analysis = null;
        analysis = Args.parseSeqAnalysis(args,"analysis");                                                        
        pvalthresh = Args.parseDouble(args,"thresh",.01);
    }
    public void getRegions(String args[]) {}
    public Collection<SeqAnalysisResult> getOverlappingRegions (Region wholeRegion) throws SQLException {
        ArrayList<SeqAnalysisResult> output = new ArrayList<SeqAnalysisResult>();
        for (SeqAnalysisResult r : analysis.getResults(getGenome(), wholeRegion)) {
            if (!Double.isInfinite(r.pvalue) && r.pvalue > pvalthresh) {
                continue;
            }
            output.add(r);
        }
        return output;
    }

}