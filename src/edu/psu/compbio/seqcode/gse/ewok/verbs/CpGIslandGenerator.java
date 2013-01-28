package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class CpGIslandGenerator<X extends Region> implements Expander<X,CpGIsland> {

    private Genome genome;
    private String tablename = "cpgIslandExt";
    
    public CpGIslandGenerator(Genome g) {
        genome = g;
    }

    public Iterator<CpGIsland> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select chromStart, chromEnd, cpgNum, gcNum, perCpg, perGc, obsExp from " + tablename + " where chrom = ? and " +
                                                        "((chromStart <= ? and chromEnd >= ?) or (chromStart >= ? and chromStart <= ?)) order by chromStart");
            String chr = region.getChrom();
            if (!chr.matches("^(chr|scaffold).*")) {
                chr = "chr" + chr;
            }
            ps.setString(1,chr);
            ps.setInt(2,region.getStart());
            ps.setInt(3,region.getStart());
            ps.setInt(4,region.getStart());
            ps.setInt(5,region.getEnd());
            ResultSet rs = ps.executeQuery();
            ArrayList<CpGIsland> results = new ArrayList<CpGIsland>();
            while (rs.next()) {
                CpGIsland island = new CpGIsland(genome,
                                                 region.getChrom(),
                                                 rs.getInt(1),
                                                 rs.getInt(2),
                                                 rs.getInt(3),
                                                 rs.getInt(4),
                                                 rs.getInt(5),
                                                 rs.getInt(6),
                                                 rs.getFloat(7));
                results.add(island);
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC CpG islands",ex);
        }
    }
}

