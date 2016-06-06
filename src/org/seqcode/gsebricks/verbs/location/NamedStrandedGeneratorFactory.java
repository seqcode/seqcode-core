/*
 * Created on Sep 28, 2006
 */
package org.seqcode.gsebricks.verbs.location;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedStrandedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.RegionExpanderFactory;
import org.seqcode.gsebricks.verbs.Expander;

public class NamedStrandedGeneratorFactory implements RegionExpanderFactory<NamedStrandedRegion> {
    String type;

    public NamedStrandedGeneratorFactory() {
        type = "NamedStrandedRegion";
    }
    public NamedStrandedGeneratorFactory(String t) {
        type = t;
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "NamedStrandedRegion";}
    public Expander<Region, NamedStrandedRegion> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, NamedStrandedRegion> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("NamedStrandedGenerator must have a type");
        } else {
            return new NamedStrandedGenerator(g, type);
        }
    }
}
