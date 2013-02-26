package edu.psu.compbio.seqcode.gse.seqview.paintable;

import edu.psu.compbio.seqcode.gse.seqview.model.ChipSeqDataModel;

import java.awt.Graphics2D;

public class SeqCombinedPainter extends SeqPainter {

    private SeqBasicOverlapPainter basic;
    private SeqAboveBelowStrandPainter stranded;

    public SeqCombinedPainter (ChipSeqDataModel model) {
        super(model);
        basic = new SeqBasicOverlapPainter(model);
        basic.setProperties(getProperties());
        stranded = new SeqAboveBelowStrandPainter(model);
        stranded.setProperties(getProperties());       
    }

    protected void paintOverlapping(Graphics2D g, 
                                    int x1, int y1, 
                                    int x2, int y2) {
        if (getProperties().Stranded) {
            stranded.paintOverlapping(g,x1,y1,x2,y2);
        } else {
            basic.paintOverlapping(g,x1,y1,x2,y2);
        }
    }
    
    protected void paintNonOverlapping(Graphics2D g, 
                                       int x1, int y1, 
                                       int x2, int y2) {
        if (getProperties().Stranded) {
            stranded.paintNonOverlapping(g,x1,y1,x2,y2);
        } else {
            basic.paintNonOverlapping(g,x1,y1,x2,y2);
        }
    }
}