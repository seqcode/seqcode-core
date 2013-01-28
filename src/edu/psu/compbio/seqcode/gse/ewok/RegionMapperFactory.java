package edu.psu.compbio.seqcode.gse.ewok;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;

public interface RegionMapperFactory<PRODUCT> {
    public void setType(String type);    
    public String getType();
    public String getProduct();
    public Mapper<Region,PRODUCT> getMapper(Genome g);
    public Mapper<Region,PRODUCT> getMapper(Genome g, String type);
}
