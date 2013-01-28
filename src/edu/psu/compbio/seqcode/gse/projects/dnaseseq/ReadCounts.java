package edu.psu.compbio.seqcode.gse.projects.dnaseseq;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;

public class ReadCounts {
    private int[] counts;
    private Region region;

    public ReadCounts(int[] c, Region r) {
        counts = c;
        region = r;
    }
    public int getCount(int pos) {
        if (pos >= region.getStart() && pos < region.getEnd()) {
            return counts[pos - region.getStart()];
        } else {
            return 0;
        }
    }
    public int[] getCounts() {return counts;}
}