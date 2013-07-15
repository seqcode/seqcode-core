package edu.psu.compbio.seqcode.projects.shaun.teqseq.analysis;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import cern.jet.random.Poisson;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

/** java edu.psu.compbio.seqcode.gse.shaun.teqseq.analysis.ThreePrimeUTR --species "$MM;mm9" --file TS.tunits.out 
 * --step 10 size step to analyze UTR in
 * --mincount 5 how many reads must cover a base to count as transcription
 * --pval .001 poisson pvalue for detecting transcription
 * --readsperbase .05 expected reads per base for poisson model
 * --threads 8  how many threads to start
 * --longestutr 10000 longest possible UTR; analysis only looks out this far
 * --analysiswindow 100  size of region in which to look for expression in each step
 */

public class ThreePrimeUTR {

    private String lastLine = null;
    private Genome genome;
    private RefGeneGenerator refgene;
    private BufferedReader reader;
    private int nCol;
    private String columnLabels[];

    private int step = 10;
    private int minCountForTranscription = 5;
    private double pval = .001;
    private double isExpressedPval = .01;
    private double diffExprPval = .01;
    private int maxThreads = 8;
    private int longestPossibleUTR = 10000;
    private int analysisWindowSize = 50;

    public static void main(String args[]) throws IOException, NotFoundException {
        ThreePrimeUTR tp = new ThreePrimeUTR(args);
        tp.execute();
    }

    public ThreePrimeUTR (String args[]) throws IOException, NotFoundException {
        reader = new BufferedReader(new FileReader(Args.parseString(args,"file",null)));        
        step = Args.parseInteger(args,"step",step);
        minCountForTranscription = Args.parseInteger(args,"mincount",minCountForTranscription);
        pval = Args.parseDouble(args,"pval",pval);
        genome = Args.parseGenome(args).cdr();
        maxThreads = Args.parseInteger(args,"threads",maxThreads);
        longestPossibleUTR = Args.parseInteger(args,"longestutr",longestPossibleUTR);
        analysisWindowSize = Args.parseInteger(args,"analysiswindow",analysisWindowSize);
        List<RefGeneGenerator> rgg = Args.parseGenes(args);
        if (rgg.size() > 0) {
            refgene = rgg.get(0);
        } else {
            refgene = new RefGeneGenerator(genome, "refGene");
        }
        refgene.retrieveExons(true);
        refgene.retrieveCoding(true);
        refgene.setWantAlias(false);
        System.err.println("Set refgene to get coding regions and not get aliases");
    }

    private void execute() throws IOException {
        Thread threads[] = new Thread[maxThreads];
        for (int i = 0; i < maxThreads; i++) {
            Thread t = new Thread(new ThreePrimeThread());
            threads[i] = t;
            t.start();
        }
        boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(10000);
            } catch (InterruptedException e) { }
            for (int i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }            
        }
    }

    class ThreePrimeThread implements Runnable {


        private Region regions[];
        private int n, biggestTranscript;
        private double counts[][];
        private double perBaseMean[];
        private Poisson poisson = new Poisson(1, new DRand());
        private Binomial binomial = new Binomial(100, .5, new DRand());
        private SequenceGenerator seqgen;
        private Pattern pattern = Pattern.compile("AATT?AAA");

        public ThreePrimeThread() {
            seqgen = new SequenceGenerator();
        }

        public void run() {
            try {
                while (fillArrays()) {
                    executeChrom();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        /** fills regions, counts, and n with the data for the next chromosome.
            returns true iff there was data.  returns false when there's no
            more data to read 
        */
        private boolean fillArrays() throws IOException {
            synchronized (reader) {
                ArrayList<String> lines = new ArrayList<String>();
                String pieces[];
                if (lastLine == null) {
                    lastLine = reader.readLine(); 
                    if (lastLine == null) {
                        return false;
                    }
                    pieces = lastLine.split("\\t");
                    nCol = pieces.length - 3;
                    columnLabels = new String[nCol];
                    for (int i = 2; i < pieces.length - 1; i++) {
                        columnLabels[i-2] = pieces[i];
                    }
                    lastLine = reader.readLine();
                }
                if (lastLine == null) { return false;}
                pieces = lastLine.split(":");
                String thischrom = pieces[0];
                System.err.println("Parsing " + thischrom);
                int tcl = thischrom.length();
                while (lastLine != null && thischrom.equals(lastLine.substring(0,tcl))) {
                    lines.add(lastLine);
                    lastLine = reader.readLine();
                }
                n = lines.size();
                regions = new Region[n];
                counts = new double[nCol][n];
                perBaseMean = new double[nCol];
                biggestTranscript = 100;
                for (int i = 0; i < n; i++) {
                    pieces = lines.get(i).split("\\t");
                    if (pieces.length != nCol + 3) {
                        System.err.println("p.l " + pieces.length + " but nCol " + nCol);
                        System.err.println("Bad line " + lastLine);
                    }

                    assert(pieces.length == nCol + 3);
                    regions[i] = Region.fromString(genome, pieces[0]);
                    if (regions[i].getWidth() > biggestTranscript) {
                        biggestTranscript = regions[i].getWidth();
                    }
                    for (int j = 2; j < pieces.length - 1; j++) {
                        counts[j-2][i] = Double.parseDouble(pieces[j]);
                    }
                }
            }
            return (n > 0);
        }
        /** getLowerBound returns the first region/line that
            needs to be considered for a window starting at pos.
        */
        private int getLowerBound(int pos) {
            int l = 0;
            int r = n;
            while (r - l > 10) {
                int c = (l+r)/2;
                if (regions[c].getStart() > pos) {
                    l = c;
                } else {
                    r = c;
                }
            }
            while (l > 0 && regions[l].getStart() > pos - biggestTranscript) {
                l--;
            }
            return l;        
        }
        /** getLowerBound returns the last region/line (inclusive) that
            needs to be considered for a window ending at pos.
        */
        private int getUpperBound(int pos) {
            int l = 0;
            int r = n;
            while (r - l > 10) {
                int c = (l+r)/2;
                if (regions[c].getStart() > pos) {
                    l = c;
                } else {
                    r = c;
                }
            }
            while (r < n && regions[r].getStart() < pos) {
                r++;
            }
            if (r == n) {
                r = n - 1;
            }
            return r;
        }
        /** returns the region that we should examine as 3' UTR for
            the gene
        */
        private Region getUTR(ExonicGene g) {
            Region lastExon = null;
            Iterator<Region> exons = g.getExons();
            while (exons.hasNext()) {
                Region e = exons.next();
                if (lastExon == null || 
                    ((g.getStrand() == '+' && e.getEnd() > lastExon.getEnd())) ||
                    (g.getStrand() == '-' && e.getStart() < lastExon.getStart())) {
                    lastExon = e;
                }
            }
            if (lastExon == null) { return null ;}
            if (g.getStrand() == '+') {
                return new Region(genome, g.getChrom(), g.getEnd(), lastExon.getEnd());
            } else {
                return new Region(genome, g.getChrom(), lastExon.getStart(), g.getStart());
            }
        }
        /* writes into lastTransPosStrict and lastTransPosLoose the
           best estimate for the last transcribed base of the utr.
           Value is the number of bases past the end of the CDS.
           Strict will be a shorter UTR and Loose will be a larger
           value (uses a looser threshold for detecting expression).
           Uses a poisson test on windows of size analysisWindowSize
           until it hits one that doesn't contain enough reads to be
           expressed according to a poisson test.
        */
        private void fillLastTransPosPoisson(ExonicGene gene, Region utr, int lastTransPosStrict[], int lastTransPosLoose[]) {
            boolean strand = gene.getStrand() == '+';
            int start = strand ? utr.getStart() : utr.getEnd();
            int direction = strand ? 1 : -1;
            for (int j = 0; j < nCol; j++) {
                lastTransPosStrict[j] = -1;
                lastTransPosLoose[j] = -1;
            }
            double sums[] = new double[nCol];
            int stillGoing = nCol;
            for (int pos = start; Math.abs(pos - start) < longestPossibleUTR; pos += direction * analysisWindowSize / 5) {
                if (stillGoing == 0) {break;}
                int first = getLowerBound(strand ? pos : pos - analysisWindowSize);
                int last = getUpperBound(strand ? pos + analysisWindowSize : pos);
                Region r = new Region(utr.getGenome(), 
                                      utr.getChrom(), 
                                      strand ? pos : pos - analysisWindowSize, 
                                      strand ? pos + analysisWindowSize : pos);
                for (int j = 0; j < nCol; j++) {sums[j] = 0;}
                for (int i = first; i <= last; i++) {
                    if (!regions[i].overlaps(r)) { continue; }
                    for (int j = 0; j < nCol; j++) {
                        sums[j] += counts[j][i];
                    }
                }
                int lenFromEnd = strand ? pos - start: start - pos;
                for (int j = 0; j < nCol; j++) {

                    double mean = analysisWindowSize * perBaseMean[j];
                    int x = (int)sums[j];
                    poisson.setMean(mean);
                    double poissonPValue = 1 - poisson.cdf(x) + poisson.pdf(x);
                    // System.err.println(String.format("%s in %d : mean %.2f, count %d, pval %.2e",
                    //                                  r.toString(), j, mean, x, poissonPValue));
                    if (poissonPValue > pval && lastTransPosLoose[j] == -1) {
                        lastTransPosLoose[j] = lenFromEnd;
                        stillGoing--;
                    }
                    if (poissonPValue * 10 > pval && lastTransPosStrict[j] == -1) {
                        lastTransPosStrict[j] = lenFromEnd;
                    }
                }
            }
        }
        /** fills expression[] with the counts of all transcriptional units that overlap the exons
         */
        private void fillGeneExpression(ExonicGene gene, double expression[]) {
            Iterator<Region> exons = gene.getExons();
            for (int j = 0; j < expression.length; j++) {
                expression[j] = 0;
            }

            while (exons.hasNext()) {
                Region exon = exons.next();
                int first = getLowerBound(exon.getStart());
                int last = getUpperBound(exon.getEnd());        
                for (int i = first; i <= last; i++) {
                    if (regions[i].overlaps(exon)) {
                        for (int j = 0; j < nCol; j++) {
                            expression[j] += counts[j][i];
                        }
                    }
                }
            }
        }
        /** returns true iff it looks like the the gene whose expression counts
         * are in expression[] (see fillGeneExpression) is differentially expressed
         * such that experiment one has a higher expression that experiment two
         */
        private boolean isDiffExpressed(double expression[], int one, int two, String label) {
            double factor = perBaseMean[one] / perBaseMean[two];
            poisson.setMean(expression[two] * factor);
            int oneval = (int)(expression[one] / 2);  // look for 2X differential expression
            double p = 1 - poisson.cdf(oneval) + poisson.pdf(oneval);
            return (p < diffExprPval);
        }
        /** returns true iff it looks like the gene is expressed in experiment 
         * i (see fillGeneExpression)
         */
        private boolean isExpressed(ExonicGene gene, double expression[], int i) {
            Iterator<Region> exons = gene.getExons();
            int exonlen = 0;
            while (exons.hasNext()) {
                Region exon = exons.next();
                exonlen += exon.getWidth();
            }
            poisson.setMean(exonlen * perBaseMean[i]);
            double p = 1 - poisson.cdf((int)expression[i]) + poisson.pdf((int)expression[i]);
            return (p < isExpressedPval);        
        }
        /* determine whether the putative transcription stop sites have the
            aatt?aaa motif w/in 20bp on either side
        */
        private void analyzeEnds(ExonicGene gene, Region utr, int longer, int shorter) {
            boolean strand = gene.getStrand() == '+';
            int mult = strand ? 1 : -1;
            int basepos = strand ? utr.getEnd() : utr.getStart();
            Region longerRegion = new Region(gene.getGenome(), gene.getChrom(),
                                             basepos + longer * mult - 20,
                                             basepos + longer * mult + 20);
            Region shorterRegion = new Region(gene.getGenome(), gene.getChrom(),
                                              basepos + shorter * mult - 20,
                                              basepos + shorter * mult + 20);
            String longerSeq = seqgen.execute(longerRegion);
            String shorterSeq = seqgen.execute(shorterRegion);
            Matcher lm = pattern.matcher(longerSeq);
            Matcher sm = pattern.matcher(shorterSeq);
            boolean lmm = lm.matches();
            boolean smm = sm.matches();
            System.out.println(String.format("\t%s vs %s",
                                             lmm ? "aataaa" : "",
                                             smm ? "aataaa" : ""));
        }
        /** 
         * run the analysis on one chromosome.  Must call fillArrays() first.
         */
        private void executeChrom() {
            if (n == 0) {
                return;
            }
            String chromName = regions[0].getChrom();
            System.err.println("Working on " + chromName);
            Region chrom = new Region(genome, chromName, 1, genome.getChromLength(chromName));
            for (int j = 0; j < nCol; j++) {
                perBaseMean[j] = 0;
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < nCol; j++) {
                    perBaseMean[j] += counts[j][i];
                }
            }
            for (int j = 0; j < nCol; j++) {
                perBaseMean[j] /= chrom.getWidth();
            }
            Iterator<Gene> genes;
            synchronized (refgene) {
                genes = refgene.execute(chrom);
            }
            int lastTransPosStrict[] = new int[nCol];
            int lastTransPosLoose[] = new int[nCol];
            double expression[] = new double[nCol];
            while (genes.hasNext()) {
                Gene g = genes.next();
                if (! (g instanceof ExonicGene)) { continue; }
                ExonicGene exonic = (ExonicGene)g;
                Region utr = getUTR(exonic);
                if (utr == null) {
                    continue;
                }
                if (utr.getWidth() < 20 || utr.getWidth() > 1000) {
                    continue;
                }
                fillGeneExpression(exonic, expression);
                fillLastTransPosPoisson(exonic, utr, lastTransPosStrict, lastTransPosLoose);
                for (int i = 0; i < nCol; i++) {
                    if (!isExpressed(exonic, expression, i)) { continue ;}
                    for (int j = 0; j < nCol; j++) {
                        if (i == j) { continue;  }
                        if (!isExpressed(exonic, expression, j)) { continue ;}
                        //                        if (isDiffExpressed(expression,i,j, g.toString())) { continue; }
                        if (lastTransPosStrict[i] - lastTransPosLoose[j] >= step) {
                            System.out.print(String.format("%s\t%s\t\t%s vs %s\tlen %d vs %d\texpr %.1f vs %.1f",
                                                           g.toString(),
                                                           utr.toString(),
                                                           columnLabels[i], columnLabels[j],
                                                           lastTransPosStrict[i], lastTransPosLoose[j],
                                                           expression[i],expression[j]));
                            analyzeEnds(exonic, utr, lastTransPosStrict[i], lastTransPosLoose[j]);
                        }
                    }
                }
            }
        }            
    }

}

