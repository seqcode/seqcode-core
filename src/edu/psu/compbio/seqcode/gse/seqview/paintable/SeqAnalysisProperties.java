package edu.psu.compbio.seqcode.gse.seqview.paintable;

public class SeqAnalysisProperties extends PaintableProperties {

    public void loadDefaults () {
        // don't load the track label from the defaults, since it varies by experiment.
        String origTrackLabel = TrackLabel;
        super.loadDefaults();
        TrackLabel = origTrackLabel;
    }

}