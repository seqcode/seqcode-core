package edu.psu.compbio.seqcode.gse.seqview.paintable;

import java.awt.Color;

public class SeqHistogramProperties extends ExperimentPaintableProperties {

    public Integer MaxReadCount = 50;
    public Boolean Stranded = Boolean.TRUE;
    public Boolean BinAutoUpdate = Boolean.FALSE;
    public Boolean DrawBinSize = Boolean.TRUE;
    public Boolean DrawPairedCurves = Boolean.FALSE;
    
    private Color plusColor = Color.red;
    private Color minusColor = Color.blue;
    private Color unstrandedColor = Color.gray;
    private Color mateArcColor = new Color(100,100,100,50);
    private Color splitReadArcColor = new Color(0,255,255,50);
    
    private int drawGaussianMaxWindow = 30000;
    private Integer LineWidth = 1;
    
    public Color getPlusColor(){return plusColor;}
    public Color getMinusColor(){return minusColor;}
    public Color getUnstrandedColor(){return unstrandedColor;}
    public Color getMateArcColor(){return mateArcColor;}
    public Color getSplitReadArcColor(){return splitReadArcColor;}
    public int getDrawGaussianMaxWindow(){return drawGaussianMaxWindow;}
    public int getLineWidth(){return LineWidth;}
}