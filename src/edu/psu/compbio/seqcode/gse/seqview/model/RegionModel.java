package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;

/* A RegionModel is a model that also carries an ewok Region with it.  */
   
public interface RegionModel extends Model {

    public void setRegion(Region r) throws NullPointerException;
    public Region getRegion();
}
