package org.seqcode.gse.gsebricks.verbs;

import org.seqcode.genome.location.Region;

public class DistanceFilter<X extends Region> implements Filter<X,X> {

    private int maxdistance;
    private Region base;

    public DistanceFilter(int d, Region b) {
        maxdistance = d;
        base = b;
    }
    public X execute(X other) {
        if (base.distance(other) < maxdistance) {
            return other;
        } else {
            return null;
        }
    }
}