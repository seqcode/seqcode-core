package org.seqcode.gse.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedTypedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.RegionExpanderFactory;
import org.seqcode.gse.gsebricks.verbs.Expander;

public class SGDOtherGeneratorFactory implements RegionExpanderFactory<NamedTypedRegion> {
    private String type;
    
    public SGDOtherGeneratorFactory() {
        type = "sgdOther";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "NamedTypedRegion";}
    public Expander<Region, NamedTypedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, NamedTypedRegion> getExpander(Genome g, String type) {
        return new SGDOtherGenerator(g,type);
    }
}
