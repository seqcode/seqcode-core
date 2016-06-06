package org.seqcode.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedTypedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.RegionExpanderFactory;
import org.seqcode.gsebricks.verbs.Expander;

public class RepeatMaskedGeneratorFactory implements RegionExpanderFactory<NamedTypedRegion> {
    private String type;
    
    public RepeatMaskedGeneratorFactory() {
        type = "RepeatMasked";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "NamedTypedRegion";}
    public Expander<Region, NamedTypedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, NamedTypedRegion> getExpander(Genome g, String type) {
        return new RepeatMaskedGenerator(g);
    }
}
