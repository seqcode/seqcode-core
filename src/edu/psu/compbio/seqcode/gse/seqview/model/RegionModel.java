package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;

/* A RegionModel is a model that also carries an ewok Region with it.  */
   
public interface RegionModel extends Model {

	//Set the region, typically ignore if region is the same as already set
    public void setRegion(Region r) throws NullPointerException;
    //Set the region even if same as already set
    public void resetRegion(Region r) throws NullPointerException;
    //Get the current Region
    public Region getRegion();
    //Check if the connection driving the model is still open.
    public boolean connectionOpen(); 
    //Reopen the connection
    public void reconnect();
}
