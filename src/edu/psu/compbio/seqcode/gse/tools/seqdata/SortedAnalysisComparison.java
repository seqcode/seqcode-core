package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.SQLException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;

/**
 * Sorts analysis two based on some criteria (pvalue, fold change, etc).
 * Runs down list of events from two until less than some percent
 * is contained in analysis one.
 *
 * [--overlap .9]  go down list two until overlap with list one drops to this percent
 * [--sort {pvalue|foldchange}]  sort analysis two on this field
 * [--fcone 1.0] minimum fold change for events from analysis one
 * [--fctwo 1.0] minimum fold change for events from analysis two
 * [--pvalone .001] max pvalue for events from analysis one
 * [--pvaltwo .001] max pvalue for events from analysis two
 * [--firstcheck 100]  accept this many events from two before performing
 *   the first overlap check with one.  This is important in case, eg, 
 *   the first event in two isn't in one; you don't want to stop there if
 *   the next 99 are in one.

 */

public class SortedAnalysisComparison extends CompareTwoAnalyses {

    private String sortType;
    private double minFoldChangeOne, minFoldChangeTwo;
    private double maxPvalueOne, maxPvalueTwo;
    private double overlapPercent = .9;
    private int firstCheck = 100;

    public SortedAnalysisComparison() {
        super();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException {
        super.parseArgs(args);
        sortType = Args.parseString(args,"sort","pvalue");
        minFoldChangeOne = Args.parseDouble(args,"fcone",1.0);
        minFoldChangeTwo = Args.parseDouble(args,"fctwo",1.0);
        maxPvalueOne = Args.parseDouble(args,"pvalone",.001);
        maxPvalueTwo = Args.parseDouble(args,"pvaltwo",.001);
        overlapPercent = Args.parseDouble(args,"overlap",.9);
        if (!(sortType.equals("pvalue") || sortType.equals("foldchange"))) {
            throw new RuntimeException("Invalid sort type "+ sortType);
        }

    }
    public List<SeqAnalysisResult> getResultsOne(Region region) throws SQLException {
        List<SeqAnalysisResult> output = new ArrayList<SeqAnalysisResult>();
        for (SeqAnalysisResult r : super.getResultsOne(region)) {
            if (r.getPValue() <= maxPvalueOne && r.getFoldEnrichment() >= minFoldChangeOne) {
                output.add(r);
            }
        }
        return output;
    }
    public List<SeqAnalysisResult> getResultsTwo(Region region) throws SQLException {
        List<SeqAnalysisResult> output = new ArrayList<SeqAnalysisResult>();
        for (SeqAnalysisResult r : super.getResultsTwo(region)) {
            if (r.getPValue() <= maxPvalueTwo && r.getFoldEnrichment() >= minFoldChangeTwo) {
                output.add(r);
            }
        }
        return output;
    }
    public List<SeqAnalysisResult> getOutputEvents() throws SQLException {
        List<SeqAnalysisResult> listOne = getResultsOne();
        Collections.sort(listOne);

        List<SeqAnalysisResult> listTwo = getResultsTwo();
        System.err.println("There are " + listOne.size() + " events in analysis one and " + listTwo.size() + " events in analysis two");
        if (sortType.equals("pvalue")) {
            Collections.sort(listTwo,new SeqAnalysisResultPvalueComparator());
        } else if (sortType.equals("foldchange")) {
            Collections.sort(listTwo,new SeqAnalysisResultEnrichmentComparator());
        } else {
            throw new RuntimeException("Invalid sort type " + sortType);
        }

        double overlap = 0;
        int i = 0;
        while (i < listTwo.size() &&
               (i < firstCheck || overlap/i > overlapPercent)) {
            if (containsMatch(listOne,listTwo.get(i))) {
                overlap++;
            }
            i++;
        }
        i--;
        System.err.println("Keeping through " + i);
        if (overlap/i > overlapPercent) {
            return listTwo.subList(0,i+1);
        } else {
            return new ArrayList<SeqAnalysisResult>();
        }
    }
    
    public static void main(String args[]) throws Exception {
        SortedAnalysisComparison sac = new SortedAnalysisComparison();
        sac.parseArgs(args);
        sac.printOutputEvents();
    }
    

}