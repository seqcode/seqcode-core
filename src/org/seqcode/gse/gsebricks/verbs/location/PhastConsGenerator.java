package org.seqcode.gse.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;

public class PhastConsGenerator<X extends Region> extends ScoredRegionGenerator<X> {
    public PhastConsGenerator(Genome g, String tablename) {
        super(g,tablename);
    }
    public String getColumnsSQL() { return "chromStart, chromEnd, sumData / (chromEnd - chromStart + 1)";}
}