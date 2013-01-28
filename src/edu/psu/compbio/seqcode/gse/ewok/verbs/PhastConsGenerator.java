package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class PhastConsGenerator<X extends Region> extends ScoredRegionGenerator<X> {
    public PhastConsGenerator(Genome g, String tablename) {
        super(g,tablename);
    }
    public String getColumnsSQL() { return "chromStart, chromEnd, sumData / (chromEnd - chromStart + 1)";}
}