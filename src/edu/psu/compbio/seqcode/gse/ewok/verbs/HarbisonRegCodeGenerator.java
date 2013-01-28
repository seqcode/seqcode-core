package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class HarbisonRegCodeGenerator<X extends Region> implements Expander<X,HarbisonRegCodeRegion> {
    
    private String tableName;
    private Genome genome;
    
    public HarbisonRegCodeGenerator(Genome g) {
        tableName = "transRegCode";
        genome = g;
    }
    
    public HarbisonRegCodeGenerator(Genome g, String tn) { 
        tableName = tn;
        genome = g;
    }

    public Iterator<HarbisonRegCodeRegion> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select name, chrom, chromStart, chromEnd, score, chipEvidence, consSpecies from " +
                    tableName + " where chrom = ? and " +
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
            ArrayList<HarbisonRegCodeRegion> results = new ArrayList<HarbisonRegCodeRegion>();
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                results.add(new HarbisonRegCodeRegion(region.getGenome(),
                                                      rs.getString(2),
                                                      rs.getInt(3),
                                                      rs.getInt(4),
                                                      rs.getString(1),
                                                      rs.getDouble(5),
                                                      rs.getString(6),
                                                      rs.getInt(7)));
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new edu.psu.compbio.seqcode.gse.utils.database.DatabaseException("Couldn't get UCSC RefGenes",ex);
        }
    }
}    
