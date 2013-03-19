package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.GeneFactory;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;

/**
 * A generator factory that is identical to RefGeneGeneratorFactory, 
 * but returns a different product name to allow a different painter to be used. 
 * @author mahony
 */
public class CDSGeneratorFactory implements RegionExpanderFactory<Gene>, GeneFactory {
    private String type;

    public CDSGeneratorFactory() {
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "CDS";}
    public Expander<Region, Gene> getExpander(Genome g) {
        return getExpander(g,type);
    }

    public Expander<Region, Gene> getExpander(Genome g, String type) {
        if (type == null) {
            return new RefGeneGenerator(g);
        } else {
            RefGeneGenerator gg = new RefGeneGenerator(g, type);
            return gg;
        }
    }
}
