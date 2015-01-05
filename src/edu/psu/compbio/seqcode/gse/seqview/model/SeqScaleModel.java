package edu.psu.compbio.seqcode.gse.seqview.model;

import java.util.EventObject;
import java.util.ArrayList;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.utils.*;

public class SeqScaleModel extends SeqViewModel implements RegionModel, Listener<EventObject> {
    
    private ArrayList<SeqDataModel> models;
    private double maxOverlap;
    private Region region;    

    public SeqScaleModel() {
        models = new ArrayList<SeqDataModel>();
        maxOverlap = 1;
    }

    public SeqScaleModel(SeqDataModel m) {
        models = new ArrayList<SeqDataModel>();
        models.add(m);
        maxOverlap = 1;
    }
    
    public void addModel(SeqDataModel m) {
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
    
    public double getMaxVal() {
        maxOverlap = 1;
        for(SeqDataModel dm : models) { 
            maxOverlap = Math.max(maxOverlap, dm.getTotalMaxOverlap());
        }
        return maxOverlap;
    }

    public void eventRegistered(EventObject o) {
        notifyListeners();
    }

}
