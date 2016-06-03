package org.seqcode.gse.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.ScoredRegion;
import org.seqcode.gse.gsebricks.RegionExpanderFactory;
import org.seqcode.gse.gsebricks.verbs.Expander;

public class ScoredRegionGeneratorFactory implements RegionExpanderFactory<ScoredRegion> {
    private String type;
    
    public ScoredRegionGeneratorFactory() {
        type = "ScoredRegion";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "ScoredRegion";}
    public Expander<Region, ScoredRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, ScoredRegion> getExpander(Genome g, String type) {
        return new ScoredRegionGenerator(g,type);
    }
}
