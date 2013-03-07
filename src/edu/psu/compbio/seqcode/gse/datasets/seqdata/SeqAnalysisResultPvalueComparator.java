package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.util.Comparator;

public class SeqAnalysisResultPvalueComparator implements Comparator<SeqAnalysisResult> {

    public boolean equals(Object o) {
        return o instanceof SeqAnalysisResultPvalueComparator;
    }
    public int compare(SeqAnalysisResult a, SeqAnalysisResult b) {
        return Double.compare(a.getPValue(),b.getPValue());
    }

}