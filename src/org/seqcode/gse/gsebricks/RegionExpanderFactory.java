package org.seqcode.gse.gsebricks;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.Expander;

public interface RegionExpanderFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Expander<Region,PRODUCT> getExpander(Genome g);
    public Expander<Region,PRODUCT> getExpander(Genome g, String type);
}
