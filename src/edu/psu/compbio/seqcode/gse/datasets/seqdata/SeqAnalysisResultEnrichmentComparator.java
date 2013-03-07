package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.util.Comparator;

public class SeqAnalysisResultEnrichmentComparator implements Comparator<SeqAnalysisResult> {

    public boolean equals(Object o) {
        return o instanceof SeqAnalysisResultEnrichmentComparator;
    }
    public int compare(SeqAnalysisResult a, SeqAnalysisResult b) {
        return Double.compare(b.getFoldEnrichment(), a.getFoldEnrichment());
    }

}