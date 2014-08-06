package edu.psu.compbio.seqcode.gse.gsebricks.verbs.location;

import java.util.regex.*;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Mapper;

/**
 * @author shaun
 */
public class StrandedPointParser implements Mapper<String,Point> {
    
    private static Pattern regPatt;
    
    static { 
        regPatt = Pattern.compile("(\\w+):(\\d+):([^:\\s]+)");
    }
    
    private Genome genome;
    private int chromIndex, startIndex, nameIndex;

    public StrandedPointParser(Genome g) {
        genome = g;
        chromIndex = 0;
        startIndex = 1;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.gsebricks.verbs.Filter#execute(null)
     */
    public StrandedPoint execute(String input) {
        String[] array = input.split("\\s+");
        String chrom = array[chromIndex];
        
        Matcher m = regPatt.matcher(chrom);
        if(m.matches()) { 
            chrom = m.group(1);
            chrom = chrom.replaceFirst("chr", "");
            int start = Integer.parseInt(m.group(2));
            char strand = '?';
            String strandstr = m.group(3);
			if(strandstr.length() > 0) { strand = strandstr.charAt(0); }
            return new StrandedPoint(genome, chrom, start, strand);
        } else { 
            System.err.println("Line \"" + input + "\" is incorrectly formatted for a StrandedRegion");
            return null;
        }
    }

}
