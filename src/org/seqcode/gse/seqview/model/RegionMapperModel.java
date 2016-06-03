package org.seqcode.gse.seqview.model;

import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.Mapper;

public class RegionMapperModel<OUT> extends MapperModel<Region,OUT> implements RegionModel {
    private Region region;

    public RegionMapperModel(Mapper<Region,OUT> ex) {
        super(ex);
    }

    public void setRegion(Region r) throws NullPointerException {
        if (r == null) {throw new NullPointerException("Region can't be null");}
        region = r;
        setInput(r);
    }
    public void resetRegion(Region r) throws NullPointerException {
        if (r == null) {throw new NullPointerException("Region can't be null");}
        region = r;
        setInput(r);
    }
    public Region getRegion() {
        return region;
    }
}
