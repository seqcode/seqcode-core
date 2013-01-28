package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class FirstEFGenerator extends RefGeneGenerator {

    public FirstEFGenerator(Genome g) {
        super(g,"firstEF");
        retrieveExons(false);
    }
    public FirstEFGenerator(Genome g, String t) {
        super(g,t);
        retrieveExons(false);
    }
    public Iterator<Gene> execute(Region region) {        
        try {            
            java.sql.Connection cxn =
                getGenome().getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select name, chrom, strand, chromStart, chromEnd " + 
                                                        " from " + getTable() + " where chrom = ? and " +
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
            Iterator<Gene> results = parseResults(ps);
            DatabaseFactory.freeConnection(cxn);
            return results;
        } catch (SQLException ex) {
            throw new edu.psu.compbio.seqcode.gse.utils.database.DatabaseException("Couldn't get UCSC RefGenes",ex);
        }
    }
}
