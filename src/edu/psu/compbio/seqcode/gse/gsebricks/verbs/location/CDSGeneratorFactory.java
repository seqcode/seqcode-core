package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.GeneFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

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
