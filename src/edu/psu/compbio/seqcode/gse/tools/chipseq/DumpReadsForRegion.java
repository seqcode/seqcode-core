package edu.psu.compbio.seqcode.gse.tools.chipseq;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.*;
import edu.psu.compbio.seqcode.gse.utils.strings.StringUtils;

/**
 * Generates a FASTA file with reads from the specified alignment in the one or more --regions 
 * specified
 *
 * java edu.psu.compbio.seqcode.gse.tools.chipseq.DumpReadsForRegion --species "$MM;mm8" --expt "exptname;exptreplicate;alignment" --region "1:100-500"
 */

public class DumpReadsForRegion {


    public static void main(String args[]) throws SQLException, NotFoundException, IOException {       
        ChipSeqLoader loader = new ChipSeqLoader();
        List<ChipSeqLocator> locators = Args.parseChipSeq(args,"expt");
        Genome genome = Args.parseGenome(args).cdr();
        List<String> alignments = new ArrayList<String>();
        for (ChipSeqLocator l : locators) {
            for (ChipSeqAlignment a : loader.loadAlignments(l,genome)) {
                alignments.add(Integer.toString(a.getDBID()));
            }
        }
        loader.close();
        List<Region> regions = Args.parseRegions(args);
        Connection cxn = DatabaseFactory.getConnection("chipseq");
        PreparedStatement stmt = cxn.prepareStatement("select id, name, sequence from chipseqreads where id in (select distinct(read) from chipseqhits " +
                                                      "where alignment in (" + StringUtils.join(alignments,",") + ") and chromosome = ? and startpos >= ? and startpos <= ? and stoppos <= ?)");

        Set<Long> seen = new HashSet<Long>();
        for (Region r : regions) {
            int chromid = r.getGenome().getChromID(r.getChrom());
            stmt.setInt(1, chromid);
            stmt.setInt(2, r.getStart());
            stmt.setInt(3, r.getEnd());
            stmt.setInt(4, r.getEnd());
            ResultSet rs = stmt.executeQuery();
            rs.setFetchSize(10000);
            while (rs.next()) {
                if (seen.contains(rs.getLong(1))) {
                    continue;
                }
                seen.add(rs.getLong(1));
                System.out.println(">" + rs.getString(2) + "\n" + rs.getString(3));
            }
            rs.close();
        }
    }
}
