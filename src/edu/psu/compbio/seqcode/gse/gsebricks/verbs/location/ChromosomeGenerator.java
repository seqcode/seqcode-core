package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import java.util.*;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

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
