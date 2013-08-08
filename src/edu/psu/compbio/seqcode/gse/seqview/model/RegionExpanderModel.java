package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;

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
