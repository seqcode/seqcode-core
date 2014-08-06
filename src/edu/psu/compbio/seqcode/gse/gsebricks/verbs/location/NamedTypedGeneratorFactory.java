/*
 * Created on Sep 28, 2006
 */
package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.NamedTypedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

public class NamedTypedGeneratorFactory implements RegionExpanderFactory<NamedTypedRegion> {
    String type;

    public NamedTypedGeneratorFactory() {
        type = "NamedTypedRegion";
    }
    public NamedTypedGeneratorFactory(String t) {
        type = t;
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "NamedTypedRegion";}
    public Expander<Region, NamedTypedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, NamedTypedRegion> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("NamedTypedGenerator must have a type");
        } else {
            return new NamedTypedGenerator(g, type);
        }
    }
}
