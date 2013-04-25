package edu.psu.compbio.seqcode.gse.tools.rnaseq;

import java.util.*;
import java.io.IOException;
import java.sql.SQLException;
import cern.jet.random.Gamma;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

/**
 * Print read counts for genes.
 *
 * Usage:
 * java DumpExpression --species "$SC;Sigmav6" --genes s288cMapped --align "Sc TAP-tagged Khd1 against 10560-6B RNA;trial 10;Novoalign trial 10 uploaded 310109"
 *
 * Output format is <gene name>\t<number of reads>
 * The number of reads reported is only reads on the same strand
 * as the gene's ORF.
 */

public class DumpExpression {

    SeqDataLoader loader;
    List<SeqAlignment> aligns;
    RefGeneGenerator genes;
    Genome genome;

    public static void main(String args[]) throws Exception {
        DumpExpression d = new DumpExpression();
        d.parseArgs(args);
        d.run();
        d.loader.close();
    }

    public DumpExpression() throws SQLException, IOException {
        loader = new SeqDataLoader();
        aligns = new ArrayList<SeqAlignment>();
    }
    public void parseArgs(String args[]) throws SQLException, NotFoundException {
        genome = Args.parseGenome(args).cdr();
        List<SeqLocator> locators = Args.parseSeqExpt(args,"align");
        for (SeqLocator locator : locators) {
            aligns.addAll(loader.loadAlignments(locator,genome));
        }
        // parseGenes returns a list of genes; just take the first one
        genes = Args.parseGenes(args).get(0);
        if (aligns.size() == 0) {
            throw new NotFoundException("--aligns didn't match any alignments");
        }
    }

    public void run() throws SQLException, IOException {
        ChromRegionIterator chroms = new ChromRegionIterator(genome);

        while (chroms.hasNext()) {
            Region chrom = chroms.next();
            Iterator<Gene> geneiter = genes.execute(chrom);
            while (geneiter.hasNext()) {
                Gene g = geneiter.next();
                int count = loader.countByRegion(aligns,g);
                System.out.println(g.toString() + "\t" + count);
            }

        }

    }
}