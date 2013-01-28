/*
 * Created on Sep 28, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.GeneFactory;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;

public class FirstEFGeneratorFactory implements RegionExpanderFactory<Gene>, GeneFactory {
    String type;

    public FirstEFGeneratorFactory() {
        type = "Gene";
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "Gene";}
    public Expander<Region, Gene> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, Gene> getExpander(Genome g, String type) {
        if (type == null) {
            return new FirstEFGenerator(g);
        } else {
            return new FirstEFGenerator(g, type);
        }
    }
}
