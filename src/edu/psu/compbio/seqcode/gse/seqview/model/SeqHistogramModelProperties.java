package edu.psu.compbio.seqcode.gse.seqview.model;


public class SeqHistogramModelProperties extends ModelProperties {

    public Integer BinWidth = 1;
    public Integer DeDuplicate = 0;
    public Boolean UseWeights = Boolean.TRUE;
    public Integer GaussianKernelVariance = 0;
    public Boolean ReadExtension = Boolean.FALSE;
    public Boolean ShowPairedReads = Boolean.FALSE;
    public Boolean ShowSingleReads = Boolean.TRUE;
}