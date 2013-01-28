package edu.psu.compbio.seqcode.gse.ewok.verbs;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;

/**
 * Maps a StrandedRegion to its five primer end.  See also <code>StrandedToThreePrime</code>
 */
public class StrandedToFivePrime<X extends StrandedRegion> implements Mapper<X,Region> {

    private int upstream, downstream;
    public StrandedToFivePrime(int up, int down) {
		upstream = up;
		downstream = down;
    }
    public Region execute(X a) {
        int start, stop;
        switch(a.getStrand()) { 
        case '+':
            start = a.getStart() - upstream;
            stop = a.getStart() + downstream;
            return new Region(a.getGenome(), a.getChrom(), start, stop);
        case '-':
            start = a.getEnd() - downstream;
            stop = a.getEnd() + upstream;
            return new Region(a.getGenome(), a.getChrom(), start, stop);
        default:
            return a;
        }
    }

}
