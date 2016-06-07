package org.seqcode.projects.seqview.model;

import java.util.*;

import org.seqcode.data.seqdata.*;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.*;


public class SeqAnalysisModel extends SeqViewModel implements RegionModel, Runnable {

    private SeqAnalysis analysis;
    private SeqAnalysisProperties props;

    private Collection<SeqAnalysisResult> results;
    private Region region;
    private boolean newinput;

    public SeqAnalysisModel(SeqAnalysis a) {
        analysis = a;
        props = new SeqAnalysisProperties();
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
    public void resetRegion(Region r) {
        if (newinput == false) {
        	region = r;
        	newinput = true;
        }
    }
    public boolean isReady() {return !newinput;}
    public Collection<SeqAnalysisResult> getResults() { return results;}
    
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
                    results = new ArrayList<SeqAnalysisResult>();
                }
                newinput = false;
                notifyListeners();

            }
        }
        
    }

}