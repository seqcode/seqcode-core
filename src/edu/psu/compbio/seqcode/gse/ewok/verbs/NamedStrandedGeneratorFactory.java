/*
 * Created on Sep 28, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedStrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;

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
