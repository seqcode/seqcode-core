package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

public class StrandedGeneratorFactory implements RegionExpanderFactory<StrandedRegion> {
    String type;

    public StrandedGeneratorFactory() {
        type = "StrandedRegion";
    }
    public StrandedGeneratorFactory(String t) {
        type = t;
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "StrandedRegion";}
    public Expander<Region, StrandedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, StrandedRegion> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("StrandedGenerator must have a type");
        } else {
            return new StrandedRegionGenerator(g, type);
        }
    }
}
