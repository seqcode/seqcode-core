/*
 * Created on Sep 28, 2006
 */
package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.NamedStrandedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

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
