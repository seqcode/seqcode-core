package edu.psu.compbio.seqcode.gse.gsebricks;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Mapper;

public interface RegionMapperFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Mapper<Region,PRODUCT> getMapper(Genome g);
    public Mapper<Region,PRODUCT> getMapper(Genome g, String type);
}
