package edu.psu.compbio.seqcode.gse.seqview.model;

public class SeqPairedEndModelProperties extends ModelProperties {

    public Double MinimumDistance=1.0;
    public Boolean DeDuplicateByPosition = false;
    public Boolean LeftAlwaysLesser = true;
    public Boolean Cluster = false;
    public Double MaxClusterDistance = -1d;

}