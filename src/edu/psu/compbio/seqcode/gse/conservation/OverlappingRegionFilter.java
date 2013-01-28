/*
 * Created on Oct 25, 2006
 */
package edu.psu.compbio.seqcode.gse.conservation;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Filter;

/**
 * @author tdanford
 */
public class OverlappingRegionFilter implements Filter<Region,Region> {
    
    private Region target;

    public OverlappingRegionFilter(Region t) {
        target = t;
    }
    
    public void setTarget(Region t) { target = t; }

    public Region execute(Region a) {
        return (target.overlaps(a)) ? a : null;
    }

}
