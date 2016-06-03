package org.seqcode.gse.gsebricks.verbs;

import org.biojava.bio.seq.*;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedRegion;


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
