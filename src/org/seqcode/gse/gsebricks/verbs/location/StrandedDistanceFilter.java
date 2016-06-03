package org.seqcode.gse.gsebricks.verbs.location;

import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gse.gsebricks.verbs.Filter;

/**
 * StrandedDistanceFilter allows through Regions that overlap the region
 * that is <code>up</code> bases upstream and </code>down</code> bases downstream of the 
 * StrandedRegion provided in the constructor.
 */
public class StrandedDistanceFilter<X extends Region> implements Filter<X,X> {

    private StrandedRegion base;

    public StrandedDistanceFilter(int up, int down, StrandedRegion b) {
        base = b.expand(up,down);
    }
    public X execute(X other) {
        if (base.overlaps(other)) {
            return other;
        } else {
            return null;
        }
    }
}