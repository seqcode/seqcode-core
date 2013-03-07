package edu.psu.compbio.seqcode.gse.projects.dnaseseq;

import java.io.*;
import java.util.*;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.motifs.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.projects.readdb.ClientException;
import edu.psu.compbio.seqcode.gse.tools.motifs.WeightMatrixScanner;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;


/**
 * "Train" an HMM without access to ChIPSeq data for the transcription factor.  Instead,
 * analyze the DNAseq data and the motif to make a guess at the HMM parameters
 *
 * java edu.psu.compbio.seqcode.gse.projects.dnaseq.FakeHMMTrain --species "$HS;hg19" --wm "CTCF;JASPAR 11/09 MA0139.1" --dnaseq "Crawford GM12878 DNaseSeq GM12878 against Input;statistical 1/11/11" --region "7:17m-18m" --cutoff .7 --minfold 1.5
 *
 * --region 1:3-5 (can give multiple)
 * --cutoff .7  motif cutoff as multiple of maximum LL score
 * --minfold 1.5 minimum fold change for hypersensitive region
 * --modelfile hmm.model  save hmm params to this file
 *
 * The analysis specified by --dnaseq isn't actually used; the
 * code just uses the same input alignments as that analysis.
 */

public class FakeHMMTrain extends HMMTrain {

    private double motifSampleFraction;
    public FakeHMMTrain() throws IOException, ClientException, SQLException {
        super();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {
        super.parseArgs(args);
        motifSampleFraction = Args.parseDouble(args,"sample",.3);
    }

    protected boolean hasBinding(int position) {
        return (Math.random() < motifSampleFraction);
    }
    public static void main(String args[]) throws Exception {
        FakeHMMTrain train = new FakeHMMTrain();
        train.parseArgs(args);
        train.train();
        train.printModel();
        train.saveModel();
    }

}