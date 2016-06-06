package org.seqcode.projects.seqview.paintable;

import org.seqcode.projects.seqview.SeqViewProperties;

public class PaintableProperties extends SeqViewProperties {

    public Boolean DrawTrackLabel = Boolean.TRUE;
    public String TrackLabel = "";
    
    public PaintableProperties() {
        super();
        DrawTrackLabel = Boolean.TRUE;
        TrackLabel = "";
    }

    public PaintableProperties (boolean drawLabel) {
        super();
        DrawTrackLabel = drawLabel;
        TrackLabel = "";
    }

    public void loadDefaults () {
        // don't load the track label from the defaults, since it varies by experiment.
        String origTrackLabel = TrackLabel;
        super.loadDefaults();
        TrackLabel = origTrackLabel;
    }
    
    public String fileSuffix() {return "svpp";}
}