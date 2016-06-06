package org.seqcode.motifs;

import java.io.*;
import java.util.*;
import java.sql.*;
import java.text.ParseException;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.DatabaseFactory;
import org.seqcode.data.connections.UnknownRoleException;
import org.seqcode.data.motifdb.*;
import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.*;
import org.seqcode.utils.*;



/**
 * Scans sequence to determine the overal density of motifs (hits / window) in windows of
 * some size.
 *
 * --cutoff .7 minimum percent (specify between 0 and 1) match to maximum motif score that counts as a match.
 * --species "$SC;SGDv2"
 * --window 100  specifies the window size to scan.  The windows are non-overlapping
 *
 * can also specify --accept to give a regex which the motif name must match
 * of --reject to specify a regex that the motif name must not match.  Remember the regex must match
 * the *entire* name, so use something like Hnf.*

 *
 * Output format is
 * chrom\tstart\tstop\tcount
 * 
 */

public class MotifDensity {
    
    public static void main(String args[]) throws Exception {
        Collection<String> accept, reject;
        Collection<WeightMatrix> matrices = new ArrayList<WeightMatrix>();
        matrices.addAll(WeightMatrix.getAllWeightMatrices());
        double cutoffpercent = Args.parseDouble(args,"cutoff",.7);
        accept = Args.parseStrings(args,"accept");
        reject = Args.parseStrings(args,"reject");
        int windowsize = Args.parseInteger(args,"window",100);
        Genome genome = Args.parseGenome(args).cdr();
        
        SequenceGenerator seqgen = new SequenceGenerator(genome);
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              genome.getDBID());
        boolean multicount = Args.parseFlags(args).contains("multicount");
        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }
                
        matrices = Args.parseWeightMatrices(args);
        System.err.println("Scanning for " + matrices.size() + " matrices");
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }
        
        List<Region> regions = Args.parseRegionsOrDefault(args);
        for (Region region : regions) {
            char[] aschars = seqgen.execute(region).toCharArray();
            int counts[] = new int[region.getWidth() / windowsize];
            int toadd[] = multicount ? null : new int[region.getWidth() / windowsize];
            for (int i = 0; i < counts.length; i++) {
                counts[i] = 0;
            }
            for (WeightMatrix m : matrices) {
                if (toadd != null) {
                    for (int i = 0; i < counts.length; i++) {
                        toadd[i] = 0;
                    }                   
                }
                List<WMHit> hits = WeightMatrixScanner.scanSequence(m,
                                                                    (float)(m.getMaxScore() * cutoffpercent),
                                                                    aschars);
                if (multicount) {
                    for (WMHit h : hits) {
                        counts[h.getStart() / windowsize]++;
                    }
                } else {
                    for (WMHit h : hits) {
                        toadd[h.getStart() / windowsize]++;
                    }
                }
                if (!multicount) {
                    for (int i = 0; i < counts.length; i++) {
                        if (toadd[i] > 0) {
                            counts[i]++;
                        }
                    }
                }
            }
            for (int i = 0; i < counts.length; i++) {
                if (counts[i] > 0) {
                    System.out.println(String.format("%s\t%d\t%d\t%d",
                                                     region.getChrom(),
                                                     region.getStart() + i * windowsize,
                                                     region.getStart() + (i+1) * windowsize,
                                                     counts[i]));
                }
            }
        }
    }
}