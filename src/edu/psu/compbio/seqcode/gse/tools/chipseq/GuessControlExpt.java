package edu.psu.compbio.seqcode.gse.tools.chipseq;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.alignments.*;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
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

    private static ChipSeqLoader loader;
    private static Collection<String> controlFactorNames;
    private static Collection<Integer> controlFactorIDs;
    private static MetadataLoader core;
    private static Genome genome;
    public static void main(String args[]) throws Exception {
        genome = Args.parseGenome(args).cdr();
        List<ChipSeqLocator> locators = Args.parseChipSeq(args);
        loader = new ChipSeqLoader(false);
        core = new MetadataLoader();
        controlFactorNames = new ArrayList<String>();
        controlFactorNames.add("WCE");
        controlFactorNames.add("IgG");
        controlFactorNames.add("Input");
        controlFactorNames.add("GFP");        
        controlFactorNames.add("null(V5)");        
        controlFactorNames.add("Control");        
        controlFactorIDs = new ArrayList<Integer>();

        boolean strict = Args.parseFlags(args).contains("strict");

        for (String f : controlFactorNames) {
            ExptTarget factor = core.findExptTarget(f);
            if (factor != null) {
                controlFactorIDs.add(factor.getDBID());
            }
        }


        Set<String> output = new HashSet<String>();
        for (ChipSeqLocator l : locators) {
            Collection<ChipSeqAlignment> alignments = new ArrayList<ChipSeqAlignment>();
            if (l.getReplicates().size() == 0) {
                alignments.addAll(loader.loadAlignments(l.getExptName(),
                                                        null,
                                                        l.getAlignName(),
                                                        null,null,null, 
                                                        genome));
            } else {
                alignments.addAll(loader.loadAlignments(l, genome));
            }
            for (ChipSeqAlignment a : alignments) {
                Collection<ChipSeqAlignment> controls = controlsForAlignment(a);
                if (controls.size() == 0) {
                    System.err.println("No controls for " + a);
                    if (strict) {
                        System.exit(1);
                    }
                }
                for (ChipSeqAlignment c : controls) {
                    output.add(c.getExpt().getName() + ";" + c.getName());
                }
            }
        }
        for (String a : output) {
            System.out.println(a);
        }

    }
    public static Collection<ChipSeqAlignment> controlsForAlignment(ChipSeqAlignment a) throws Exception {
        ArrayList<ChipSeqAlignment> output = new ArrayList<ChipSeqAlignment>();
        ChipSeqExpt expt = a.getExpt();
        CellLine cells = expt.getCells();
        ExptCondition cond = expt.getCondition();
        for (Integer i : controlFactorIDs) {
            output.addAll(loader.loadAlignments(null,null,null,
                                                i, cells.getDBID(), cond.getDBID(), 
                                                genome));
        }
        return output;
    }
    

}