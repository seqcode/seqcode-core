package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.RepeatMaskedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 * Maps a Region to RepeatMaskedRegions from a UCSC annotations table.
 */
public class RepeatMaskedGenerator<X extends Region> implements Expander<X,RepeatMaskedRegion> {

    private Genome genome;
    private boolean onetable;
    private static Map<String,String> genomeToTablesuffix;

    static { 
        genomeToTablesuffix = new HashMap<String,String>();
        ResourceBundle res = ResourceBundle.getBundle("edu.psu.compbio.seqcode.gse.ewok.rmsk");
        Enumeration<String> keys = res.getKeys();
        while(keys.hasMoreElements()) { 
            String key = keys.nextElement();
            String value = res.getString(key);
            genomeToTablesuffix.put(key, value);
        }
    }

    public RepeatMaskedGenerator(Genome g) {
        genome = g;
        onetable = false;
    }

    public Iterator<RepeatMaskedRegion> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            String chr = region.getChrom();
            if (!chr.matches("^(chr|scaffold).*")) {
                chr = "chr" + chr;
            }
            String tablesuffix = "rmsk";
            if(genomeToTablesuffix.containsKey(genome.getVersion()))
                tablesuffix = genomeToTablesuffix.get(genome.getVersion());
                
            String tablename = onetable ? tablesuffix : (chr + "_" + tablesuffix);
            PreparedStatement ps = cxn.prepareStatement("select genoStart, genoEnd, strand, repName, repClass, repFamily, swScore from " + tablename + 
                                                        " where genoName = ? and ((genoStart <= ? and genoEnd >= ?) or (genoStart >= ? and genoStart <= ?)) order by genoStart");
            ps.setString(1,chr);
            ps.setInt(2,region.getStart());
            ps.setInt(3,region.getStart());
            ps.setInt(4,region.getStart());
            ps.setInt(5,region.getEnd());
            ResultSet rs = ps.executeQuery();
            ArrayList<RepeatMaskedRegion> results = new ArrayList<RepeatMaskedRegion>();
            while (rs.next()) {
                RepeatMaskedRegion r = new RepeatMaskedRegion(genome,region.getChrom(),
                                                              rs.getInt(1),
                                                              rs.getInt(2),
                                                              rs.getString(4),
                                                              rs.getString(5),
                                                              rs.getString(6),
                                                              (double)rs.getInt(7),
                                                              rs.getString(3).charAt(0));
                results.add(r);
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            if (!onetable) {
                onetable = true;
                return execute(region);
            } else {
                ex.printStackTrace();
                throw new DatabaseException("Couldn't get UCSC repeat masked regions",ex);
            }
        }
    }
}

