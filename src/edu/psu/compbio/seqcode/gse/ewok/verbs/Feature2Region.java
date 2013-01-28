package edu.psu.compbio.seqcode.gse.ewok.verbs;

import org.biojava.bio.seq.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;

public class Feature2Region implements Mapper<Feature,NamedRegion> {
    private Genome g;
    
    public Feature2Region(Genome g) {
        this.g = g;
    }
    
    public NamedRegion execute (Feature f) {
        return new NamedRegion(g,
                          f.getSequence().getName(),
                          f.getLocation().getMin(),
                          f.getLocation().getMax(), f.getType());
    }
}
