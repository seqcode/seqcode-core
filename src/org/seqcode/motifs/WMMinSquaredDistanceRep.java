package org.seqcode.motifs;

import java.util.*;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.ml.clustering.Cluster;
import org.seqcode.ml.clustering.ClusterRepresentative;


public class WMMinSquaredDistanceRep implements ClusterRepresentative<WeightMatrix> {
    
    private WMComparator comp;
    public WMMinSquaredDistanceRep (WMComparator c) {
        comp = c;
    }
    public WeightMatrix getRepresentative(Cluster<WeightMatrix> cluster) {
        Set<WeightMatrix> matrices = cluster.getElements();        
        WeightMatrix bestwm = null;
        double bestdist = Double.MAX_VALUE;
        for (WeightMatrix i : matrices) {
            double sum = 0;
            for (WeightMatrix j : matrices) {
                sum += Math.pow(comp.compare(i,j),2);
            }
            //            System.err.println("  " + i + " : " + sum + " <? " + bestdist);
            sum = sum / matrices.size();
            if (sum < bestdist) {
                bestwm = i;
                bestdist = sum;
            }
        }
        if (bestwm == null) {
            System.err.println("OOPS!"  + bestdist);
            System.err.println(matrices.toString());
        }
        return bestwm;
    }

}
