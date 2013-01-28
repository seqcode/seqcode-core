package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;
import java.util.regex.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;

public class PointParser implements Mapper<String,Point> {
    
    private static Pattern pPatt;
    
    static { 
        pPatt = Pattern.compile("(\\w+):(\\d+)");
    }
    
    private Genome genome;
    private int chromIndex, startIndex, minLength;

    public PointParser(Genome g) {
        genome = g;
        chromIndex = 0;
        startIndex = 1;
        minLength = startIndex + 1;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(null)
     */
    public Point execute(String input) {
        String[] array = input.split("\\s+");
        String chrom = array[chromIndex];
        chrom = chrom.replaceFirst("chr", "");
        
        Matcher m = pPatt.matcher(chrom);
        if(m.matches()) { 
            chrom = m.group(1);
            int start = Integer.parseInt(m.group(2));
            return new Point(genome, chrom, start);
        } else { 
            if(array.length >= minLength) { 
                int start = Integer.parseInt(array[startIndex]);
                return new Point(genome, chrom, start);
            } else { 
                System.err.println("Line \"" + input + "\" doesn't have the correct length (" + minLength + ")");
                return null;
            }
        }
    }
}
