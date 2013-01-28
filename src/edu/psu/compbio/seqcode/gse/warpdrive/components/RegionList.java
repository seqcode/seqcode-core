package edu.psu.compbio.seqcode.gse.warpdrive.components;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;

public interface RegionList {
    public void addRegion(Region r);
    public int regionListSize();
    public Region regionAt(int i);
}
