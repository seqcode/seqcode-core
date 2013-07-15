package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;

/* gives an iterator over all of the chromosomes in a genome */

public class ChromosomeGenerator<X extends Genome> implements Expander<X,Region> {
    
    public Iterator<Region> execute (X genome) {
        List<String> names = genome.getChromList();
        List<Region> chroms = new ArrayList<Region>();
        for (int i = 0; i < names.size(); i++) {
            chroms.add(new Region(genome,
                                  names.get(i),
                                  1,
                                  genome.getChromLength(names.get(i))));
        }
        return chroms.iterator();
    }
}
