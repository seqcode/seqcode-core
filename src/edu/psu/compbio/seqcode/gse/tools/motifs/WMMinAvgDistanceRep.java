package edu.psu.compbio.seqcode.gse.tools.motifs;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.clustering.*;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.*;

public class WMMinAvgDistanceRep implements ClusterRepresentative<WeightMatrix> {
    
    private WMComparator comp;
    public WMMinAvgDistanceRep (WMComparator c) {
        comp = c;
    }
    public WeightMatrix getRepresentative(Cluster<WeightMatrix> cluster) {
        Set<WeightMatrix> matrices = cluster.getElements();        
        WeightMatrix bestwm = null;
        double bestdist = Double.MAX_VALUE;
        for (WeightMatrix i : matrices) {
            double sum = 0;
            for (WeightMatrix j : matrices) {
                sum += comp.compare(i,j);
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
