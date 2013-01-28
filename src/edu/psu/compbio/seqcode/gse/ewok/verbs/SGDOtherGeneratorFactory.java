package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedTypedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;

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
