package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

public class RegionExpanderModel<OUT> extends ExpanderModel<Region,OUT> implements RegionModel {
    private Region region;

    public RegionExpanderModel(Expander<Region,OUT> ex) {
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
    public boolean connectionOpen(){return true;}
    public void reconnect(){}

}
