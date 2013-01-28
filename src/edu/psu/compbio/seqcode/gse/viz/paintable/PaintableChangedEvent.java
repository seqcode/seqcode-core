package edu.psu.compbio.seqcode.gse.viz.paintable;

public class PaintableChangedEvent {
    
    private Paintable fSource;
    public PaintableChangedEvent(Paintable p) { 
        fSource = p; 
    }
    
    public Paintable getSource() { return fSource; }
}

