package edu.psu.compbio.seqcode.gse.ewok;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;

public interface RegionExpanderFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Expander<Region,PRODUCT> getExpander(Genome g);
    public Expander<Region,PRODUCT> getExpander(Genome g, String type);
}
