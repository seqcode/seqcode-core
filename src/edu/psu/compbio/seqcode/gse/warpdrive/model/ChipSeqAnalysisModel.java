package edu.psu.compbio.seqcode.gse.warpdrive.model;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.*;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.utils.*;

public class ChipSeqAnalysisModel extends WarpModel implements RegionModel, Runnable {

    private ChipSeqAnalysis analysis;
    private ChipSeqAnalysisProperties props;

    private Collection<ChipSeqAnalysisResult> results;
    private Region region;
    private boolean newinput;

    public ChipSeqAnalysisModel(ChipSeqAnalysis a) {
        analysis = a;
        props = new ChipSeqAnalysisProperties();
    }
    public void clearValues() {
        results = null;
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public boolean isReady() {return !newinput;}
    public Collection<ChipSeqAnalysisResult> getResults() { return results;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    results = analysis.getResults(region);
                } catch (Exception e) {
                    e.printStackTrace();
                    results = new ArrayList<ChipSeqAnalysisResult>();
                }
                newinput = false;
                notifyListeners();

            }
        }
        
    }

}