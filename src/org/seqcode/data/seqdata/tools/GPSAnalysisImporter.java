package org.seqcode.data.seqdata.tools;

import java.io.IOException;
import java.sql.SQLException;
import java.util.*;

import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.seqdata.SeqAnalysisResult;
import org.seqcode.gsebricks.verbs.chipseq.GPSParser;
import org.seqcode.gsebricks.verbs.chipseq.GPSPeak;
import org.seqcode.utils.NotFoundException;


/**
 * See AnalysisImporter docs.  Command line options are the same; the only difference
 * is that GPSAnalysisImporter parses the GPS native output format.
 */

public class GPSAnalysisImporter extends AnalysisImporter {

    /* oracle complains about underflow if we don't limit the pvalues.  the actual 
       min value is somewhere between E-100 and E-200, but I didn't bother tracking 
       it down more closely since I don't think the difference really matters
    */
    public final static double minpval = Math.pow(10,-100);
    private Set<String> seenPositions = new HashSet<String>();

    private int lineno = 0;

    public static void main(String args[]) throws NotFoundException, SQLException, DatabaseException, IOException {
        GPSAnalysisImporter importer = new GPSAnalysisImporter();
        importer.parseArgs(args);
        importer.run(System.in);
        importer.close();
    }
    public SeqAnalysisResult parseLine(String line) {
        if (line.matches("^Position.*")) {
            return null;
        }

        GPSPeak p = GPSParser.parseLine(getGenome(),
                                        line,
                                        ++lineno);

        String k = p.getChrom() + p.getLocation();
        if (seenPositions.contains(k)) {
            return null;
        }
        seenPositions.add(k);

        return new SeqAnalysisResult(getGenome(),
                                         p.getChrom(),
                                         p.getLocation(),
                                         p.getLocation()+1,
                                         p.getLocation(),
                                         p.getStrength(),
                                         p.getControlStrength(),
                                         p.getStrength(),
                                         p.getShape(),
                                         Math.max(p.getPvalue(), minpval),
                                         p.getStrength()/p.getControlStrength());
    }


}