package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.SQLException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;

/**
 * Just report on the overlap between the two experiments
 */

public class SimpleAnalysisOverlap extends CompareTwoAnalyses {

    public static void main(String args[]) throws Exception {
        SimpleAnalysisOverlap sao = new SimpleAnalysisOverlap();
        sao.parseArgs(args);
        sao.printReport();
    }
    public void printReport() throws SQLException {
        List<SeqAnalysisResult> one = getResultsOne();
        List<SeqAnalysisResult> two = getResultsTwo();
        Collections.sort(one);
        int onefound = 0, twofound = 0;
        for (SeqAnalysisResult r : two) {
            if (containsMatch(one, r)) {
                twofound++;
            }
        }
        Collections.sort(two);
        for (SeqAnalysisResult r : one) {
            if (containsMatch(two, r)) {
                onefound++;
            }
        }
        long genomeSize = getGenome().getGenomeSize();
        double bins = genomeSize / getMaxDistance();
        double pone = one.size() / bins;
        double ptwo  = two.size() / bins;
        //        System.err.println("Pone is " + pone + " and ptwo is " + ptwo);
        double expone = one.size() * ptwo;
        double exptwo = two.size() * pone;
        System.out.println(String.format("One: %d out of %d.  Expected %.2e.  Enrichment %.1f\nTwo: %d out of %d.  Expected %.2e.  Enrichment %.1f",
                                         onefound, one.size(), expone, onefound / expone,
                                         twofound, two.size(), exptwo, twofound / exptwo));
    }
    public List<SeqAnalysisResult> getOutputEvents() { return null;}


}