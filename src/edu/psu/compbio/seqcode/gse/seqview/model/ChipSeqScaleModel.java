package edu.psu.compbio.seqcode.gse.seqview.model;

import java.util.EventObject;
import java.util.ArrayList;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.GenericExperiment;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.utils.*;

public class ChipSeqScaleModel extends SeqViewModel implements RegionModel, Listener<EventObject> {
    
    private ArrayList<ChipSeqDataModel> models;
    private double maxOverlap;
    private Region region;    

    public ChipSeqScaleModel() {
        models = new ArrayList<ChipSeqDataModel>();
        maxOverlap = 1;
    }

    public ChipSeqScaleModel(ChipSeqDataModel m) {
        models = new ArrayList<ChipSeqDataModel>();
        models.add(m);
        maxOverlap = 1;
    }
    
    public void addModel(ChipSeqDataModel m) {
        models.add(m);
        m.addEventListener(this);
    }
        
    public boolean isReady() {
        boolean isready = true;
        for (int i = 0; i < models.size(); i++) {
            isready = isready && models.get(i).isReady();
        }
        return isready;
    }
    
    public void setRegion(Region r) {
        maxOverlap = 1;
        region = r;
    }
    public void resetRegion(Region r) {
        maxOverlap = 1;
        region = r;
    }
    
    public Region getRegion() {return region;}
    public boolean connectionOpen(){return true;}
    public void reconnect(){}
    
    public double getMaxVal() {
        maxOverlap = 1;
        for(ChipSeqDataModel dm : models) { 
            maxOverlap = Math.max(maxOverlap, dm.getTotalMaxOverlap());
        }
        return maxOverlap;
    }

    public void eventRegistered(EventObject o) {
        notifyListeners();
    }

}
