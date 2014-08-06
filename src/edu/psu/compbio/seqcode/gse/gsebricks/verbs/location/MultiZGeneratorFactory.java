package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.ScoredRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

public class MultiZGeneratorFactory implements RegionExpanderFactory<ScoredRegion> {
    private String type;
    
    public MultiZGeneratorFactory() {
        type = "MultiZ";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "ScoredRegion";}
    public Expander<Region, ScoredRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, ScoredRegion> getExpander(Genome g, String type) {
        return new MultiZGenerator(g);
    }
}
