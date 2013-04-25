package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.alignments.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;



/**
 * Tries to determine one or more control experiments to use for
 * one or more chipseq experiments.
 *
 * Usage:
 * java edu.psu.compbio.seqcode.gse.tools.chipseq.GuessControlExpt --species 'Mus musculus;mm9' --chipseq 'Sing ES CTCF E14;bowtie_unique'
 */
public class GuessControlExpt {

    private static SeqDataLoader loader;
    private static Collection<String> controlTypeNames;
    private static Collection<Integer> controlTypeIDs;
    private static MetadataLoader core;
    private static Genome genome;
    public static void main(String args[]) throws Exception {
        genome = Args.parseGenome(args).cdr();
        List<SeqLocator> locators = Args.parseSeqExpt(args);
        loader = new SeqDataLoader(false);
        core = new MetadataLoader();
        controlTypeNames = new ArrayList<String>();
        controlTypeNames.add("CONTROL");
        controlTypeNames.add("INPUT");
        controlTypeIDs = new ArrayList<Integer>();

        boolean strict = Args.parseFlags(args).contains("strict");

        for (String f : controlTypeNames) {
            ExptType type = core.findExptType(f);
            if (type != null) {
                controlTypeIDs.add(type.getDBID());
            }
        }


        Set<String> output = new HashSet<String>();
        for (SeqLocator l : locators) {
            Collection<SeqAlignment> alignments = new ArrayList<SeqAlignment>();
            if (l.getReplicates().size() == 0) {
                alignments.addAll(loader.loadAlignments(l.getExptName(),
                                                        null,
                                                        l.getAlignName(),
                                                        null,null,null, 
                                                        null,null,null, 
                                                        genome));
            } else {
                alignments.addAll(loader.loadAlignments(l, genome));
            }
            for (SeqAlignment a : alignments) {
                Collection<SeqAlignment> controls = controlsForAlignment(a);
                if (controls.size() == 0) {
                    System.err.println("No controls for " + a);
                    if (strict) {
                        System.exit(1);
                    }
                }
                for (SeqAlignment c : controls) {
                    output.add(c.getExpt().getName() + ";" + c.getName());
                }
            }
        }
        for (String a : output) {
            System.out.println(a);
        }

    }
    public static Collection<SeqAlignment> controlsForAlignment(SeqAlignment a) throws Exception {
        ArrayList<SeqAlignment> output = new ArrayList<SeqAlignment>();
        SeqExpt expt = a.getExpt();
        CellLine cells = expt.getCellLine();
        ExptCondition cond = expt.getExptCondition();
        for (Integer i : controlTypeIDs) {
            output.addAll(loader.loadAlignments(null,null,null,
                                                i, null, cond.getDBID(),
                                                null, cells.getDBID(), null,
                                                genome));
        }
        return output;
    }
    

}