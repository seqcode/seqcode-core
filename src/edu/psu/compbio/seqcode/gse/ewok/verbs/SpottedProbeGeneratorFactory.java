package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.SpottedProbe;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;

public class SpottedProbeGeneratorFactory implements RegionExpanderFactory<SpottedProbe> {
    String table;

    public SpottedProbeGeneratorFactory() {
        table = "spottedArray";
    }
    public SpottedProbeGeneratorFactory(String t) {
        table = t;
    }

    public void setType(String t) {table = t;}
    public String getType() {return table;}
    public String getProduct() {return "SpottedProbe";}
    
    public Expander<Region,SpottedProbe> getExpander(Genome g) {
        return getExpander(g,table);
    }

    public Expander<Region,SpottedProbe> getExpander(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("SpottedProbeGenerator must have a type");
        } else {
            return new SpottedProbeGenerator(g, table);
        }
    }
}
