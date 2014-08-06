package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import java.util.*;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 *  Maps a region to a set of scored regions that describe the
 *  aggregate MultiZ alignment score from the UCSC annotations.
*/

public class MultiZGenerator<X extends Region> extends ScoredRegionGenerator {
    
    private static Map<String,String> genomeToTables;
    
    static { 
        genomeToTables = new HashMap<String,String>();
        ResourceBundle res = ResourceBundle.getBundle("edu.psu.compbio.seqcode.gse.gsebricks.multiz");
        Enumeration<String> keys = res.getKeys();
        while(keys.hasMoreElements()) { 
            String key = keys.nextElement();
            String value = res.getString(key);
            genomeToTables.put(key, value);
        }
    }
    
    public MultiZGenerator(Genome g) {
        super(g,null);
        String tablename;
        String version = g.getVersion();
        
        if(genomeToTables.containsKey(g.getVersion())) { 
            tablename = genomeToTables.get(g.getVersion());
        } else {
            throw new DatabaseException("Don't know the name of the multiz table for " + version);
        }
        setTablename(tablename);
    }
}

