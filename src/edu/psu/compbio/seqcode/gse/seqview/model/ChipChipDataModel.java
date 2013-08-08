package edu.psu.compbio.seqcode.gse.seqview.model;

import java.util.Hashtable;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.GenericExperiment;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.seqview.paintable.SeqViewPaintable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

/* database thread for GenericExperiment (parent of ChipChipData, ChipChipMSP, ChipChipMLE, ChipChipBayes.
   This class starts a new thread that waits on the parent.  The parent calls setRegion()
   and then notifies this thread.  This thread calls the GenericExperiment.window() method
   to start retrieving data from the database.   */

public class ChipChipDataModel extends SeqViewModel implements Runnable, RegionModel {

    private GenericExperiment data;
    private boolean newregion;
    private Region region;
    private ChipChipDataProperties props;

    public ChipChipDataModel(GenericExperiment data) {
        super();
        this.data = data;
        newregion = false;
        props = new ChipChipDataProperties();
    }
    public ChipChipDataProperties getProperties () {
        return props;
    }

    public synchronized void run() {
        while (keepRunning()) {
            try {
                if (!newregion) {
                    wait();
                }
            } catch (InterruptedException ex) {
                
            }
            if (newregion) {
                if (data instanceof SQLData) {
                    SQLData s = (SQLData) data;
                    s.regularize(getProperties().RegularizationConst);
                    s.setMaxCount(getProperties().MaxGenomicLocations);
                }
                try {
                    data.window(region.getChrom(),region.getStart(),region.getEnd());
                    newregion = false;
                    notifyListeners();                    
                } catch (NotFoundException ex) {
                    // don't do anything here.  If we got an invalid chromosome, nothing
                    // will happen and the user will notice (hopefully)
                    ex.printStackTrace();
                }
            }
        }
    }

    public synchronized void setRegion(Region r) {
        if (newregion == false) {
            newregion = true;
            region = r;
        }
    }
    public synchronized void resetRegion(Region r) {
        if (newregion == false) {
            newregion = true;
            region = r;
        }
    }
    public Region getRegion() {return region;}
    public boolean connectionOpen(){return true;}
    public void reconnect(){}
    public boolean isReady() {return !newregion;}
    public GenericExperiment getGenericExperiment() {
        return data;
    }

}

