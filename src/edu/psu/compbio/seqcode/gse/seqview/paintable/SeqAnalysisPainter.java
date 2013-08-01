package edu.psu.compbio.seqcode.gse.seqview.paintable;

import java.awt.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.seqview.model.ChipSeqAnalysisModel;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.viz.DynamicAttribute;
import edu.psu.compbio.seqcode.gse.viz.colors.ColorSet;

public class SeqAnalysisPainter extends RegionPaintable {


    private SeqAnalysis analysis;
    private ChipSeqAnalysisModel model;
    private SeqAnalysisProperties props;
    private DynamicAttribute attrib;
    private NonOverlappingLayout<SeqAnalysisResult> layout;
    private ColorSet cs;

    public SeqAnalysisPainter(SeqAnalysis a, ChipSeqAnalysisModel m) {
        super();
        analysis = a;
        model = m;
        model.addEventListener(this);
        props = new SeqAnalysisProperties();
        attrib = DynamicAttribute.getGlobalAttributes();
        layout = new NonOverlappingLayout<SeqAnalysisResult>();
        cs = new ColorSet();
        initLabels();
    }

    public SeqAnalysisProperties getProperties() {return props;}

    public int getMaxVertSpace() { 
        int numTracks = layout.getNumTracks();
        return Math.min(Math.max(40,numTracks * 12),120);
    }
    public int getMinVertSpace() { 
        int numTracks = layout.getNumTracks();
        return Math.min(Math.max(40,numTracks * 12),120);
    }
    public boolean canPaint() {
        return model.isReady();
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }
    public synchronized void eventRegistered(EventObject e) {        
        if (e.getSource() == model && model.isReady()) {
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        if (!canPaint()) {
            return;
        }
        Collection<SeqAnalysisResult> results = model.getResults();

        layout.setRegions(results);

        int numTracks = layout.getNumTracks();
        int w = lrx - ulx, h = lry - uly;
        int trackHeight = numTracks > 0 ? Math.max(1, h / numTracks) : 1;
        int spacing = Math.max(2, trackHeight/10);
        Region region = model.getRegion();
        int start = region.getStart(), end = region.getEnd();

        cs.reset();
        clearLabels();
        for (SeqAnalysisResult r : results) {
            int x1 = getXPos(r.getStart(),
                             region.getStart(),
                             region.getEnd(),
                             ulx, lrx);
            int x2 = getXPos(r.getEnd(),
                             region.getStart(),
                             region.getEnd(),
                             ulx, lrx);
            if (x2 == x1) {
                x2 = x1 + 1;
            }

            int track = layout.getTrack(r);
            int y1 = uly + trackHeight * track;            
            g.setColor(cs.getColor());
            g.fillRect(x1,y1,x2-x1,trackHeight-spacing);
            addLabel(x1,y1,x2-x1,trackHeight-spacing, String.format("%.1f/%.1f p=%.2f",
                                                                    r.foregroundReadCount,
                                                                    r.backgroundReadCount,
                                                                    r.pvalue));
        }
        if (getProperties().DrawTrackLabel) {
            g.setFont(attrib.getLargeLabelFont(w,h));
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),ulx + g.getFont().getSize()*2,uly + g.getFont().getSize());
        }
    }

}
