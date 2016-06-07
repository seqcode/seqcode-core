package org.seqcode.motifs;

import java.util.*;
import java.sql.*;

import org.seqcode.data.motifdb.*;
import org.seqcode.gseutils.*;


/**
 * For each input motif, searches entire database for similar motifs.  
 * Produces an output image for each showing the input motif and
 * those that are similar
 */

public class SimilarMotifs {

    
    public static void main(String args[]) throws NotFoundException, SQLException {
        boolean normalize = Args.parseFlags(args).contains("normalize");
        int compareLength = Args.parseInteger(args,"compareLength",-1);
        double maxDistance = Args.parseDouble(args,"maxDistance",3.0);

        WMComparator comparator = new WMDistanceComparator(normalize,compareLength);

        Collection<WeightMatrix> inputmatrices = Args.parseWeightMatrices(args);
        System.err.println("Looking for matches to " + inputmatrices.size() + " matrices");
        Collection<WeightMatrix> allmatrices = WeightMatrix.getAllWeightMatrices();
            
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              Args.parseGenome(args).cdr().getDBID());
        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }

        for (WeightMatrix m : allmatrices) {
            m.toFrequency(bgModel);
        }
        for (WeightMatrix m : inputmatrices) {
            m.toFrequency(bgModel);
        }
            
        for (WeightMatrix m : inputmatrices) {
            ArrayList<WeightMatrix> cluster = new ArrayList<WeightMatrix>();
            cluster.add(m);
            for (WeightMatrix other : allmatrices) {
                if (other.equals(m)) { continue; }
                double distance = comparator.compare(m, other);
                if (distance <= maxDistance) {
                    cluster.add(other);
                }
            }
            if (cluster.size() > 1) {
                String fname = m.toString().replaceAll("[^\\w\\d]","_") + ".png";
                ClusterMotifs.drawCluster(cluster,fname);
                for (WeightMatrix wm : cluster) {
                    System.out.println(wm.getName() + "\t" + wm.getVersion());
                }
            }
        }
    }


}