package edu.psu.compbio.seqcode.gse.ewok.verbs.chippet;

import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.chippet.ChipPetDatum;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.RegionExpanderFactory;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.Interval;

public class ChipPetExpanderFactory implements RegionExpanderFactory<ChipPetDatum> {
    
    private String name;
    
    public ChipPetExpanderFactory() {
        name = "";
    }

    public void setType(String t) {name = t;}
    
    public String getType() {return name;}
    
    public String getProduct() {return "ChipPet";}
    
    public Expander<Region,ChipPetDatum> getExpander(Genome g) {
        return getExpander(g, name);
    }

    public Expander<Region,ChipPetDatum> getExpander(Genome g, String name) {
        try {
            return new ChipPetExpander(name);
        } catch (SQLException e) {
            e.printStackTrace();
            throw new IllegalArgumentException();
        }
    }

}
