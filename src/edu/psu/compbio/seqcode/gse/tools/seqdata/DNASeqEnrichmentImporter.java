package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.IOException;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAnalysisResult;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.GPSParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.GPSPeak;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;

/**
 * See AnalysisImporter docs.  Command line options are the same; the only difference
 * is that DNASeqEnrichmentImporter parses output from DNASeqEnrichmentCaller
 */

public class DNASeqEnrichmentImporter extends AnalysisImporter {

    /* oracle complains about underflow if we don't limit the pvalues.  the actual 
       min value is somewhere between E-100 and E-200, but I didn't bother tracking 
       it down more closely since I don't think the difference really matters
    */
    public final static double minpval = Math.pow(10,-100);

    private int lineno = 0;

    public static void main(String args[]) throws NotFoundException, SQLException, DatabaseException, IOException {
        DNASeqEnrichmentImporter importer = new DNASeqEnrichmentImporter();
        importer.parseArgs(args);
        importer.run(System.in);
        importer.close();
    }
    public SeqAnalysisResult parseLine(String line) {
        if (lineno++ == 0) {
            if (!line.equals("Chrom\tStart\tEnd\tCenter\tFG\tBG\tRatio\tPvalue")) {
                throw new RuntimeException("Invalid header line: " + line);
            }
            return null;
        }
        String pieces[] = line.split("\t");
        String pospieces[] = pieces[0].split("[\\:\\-]");
        String peakpieces[] = pieces[2].split(":");
        double pval = Double.parseDouble(pieces[6]);
        if (pval < minpval) {
            pval = minpval;
        }

        return new SeqAnalysisResult(getGenome(),
                                         pieces[0],
                                         Integer.parseInt(pieces[1]),
                                         Integer.parseInt(pieces[2]),
                                         Integer.parseInt(pieces[3]),
                                         Double.parseDouble(pieces[4]),
                                         Double.parseDouble(pieces[5]),
                                         0.0,0.0,
                                         Double.parseDouble(pieces[7]),
                                         Double.parseDouble(pieces[6]));

    }


}