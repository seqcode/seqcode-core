package edu.psu.compbio.seqcode.gse.seqview.paintable;

import edu.psu.compbio.seqcode.gse.seqview.SeqViewProperties;


public class SeqDataBatchProperties extends SeqViewProperties {

	//For SeqHistogramProperties
    public Integer MaxReadCount = 50;
    public Boolean DrawBinSize = Boolean.TRUE;
    public Boolean DrawTrackLabel = Boolean.TRUE;
    
    //For SeqHistogramModelProperties
    public Integer BinWidth = 1;
    public Integer DeDuplicate = 0;
    public Boolean UseWeights = Boolean.TRUE;
    public Integer GaussianKernelVariance = 0;
    public Boolean ReadExtension = Boolean.FALSE;
    
    //For InteractionArcModelProperties
    public Integer ArcDeDuplicate = 1;

    public void loadDefaults () {
        super.loadDefaults();
    }
    
    public String fileSuffix() {return "svpp";}
}