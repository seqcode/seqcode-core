package org.seqcode.gsebricks;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.Mapper;

public interface RegionMapperFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Mapper<Region,PRODUCT> getMapper(Genome g);
    public Mapper<Region,PRODUCT> getMapper(Genome g, String type);
}
