package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.sql.SQLException;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqHit;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.Interval;

public class PointExpanderFactory implements RegionExpanderFactory<Point> {
    
    private String name;
    
    public PointExpanderFactory() {
        name = "";
    }

    public void setType(String t) {name = t;}
    
    public String getType() {return name;}
    
    public String getProduct() {return "Point";}
    
    public Expander<Region,Point> getExpander(Genome g) {
        return getExpander(g, name);
    }

    public Expander<Region,Point> getExpander(Genome g, String name) {
        try {
            if(name.startsWith("albert")) { 
                return new H2AZMidpointExpander(exptName, alignName);
            } else { 
                throw new IllegalArgumentException();
            }
        } catch (SQLException e) {
            e.printStackTrace();
            throw new IllegalArgumentException();
        }
    }
    public static String exptName = "Albert H2A.Z ChIP-Seq";
    public static String alignName = "0.90 Albert Alignment";    
}
