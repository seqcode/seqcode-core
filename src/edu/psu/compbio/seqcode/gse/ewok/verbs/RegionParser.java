/*
 * Created on Mar 9, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;
import java.util.regex.*;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;

/**
 * @author tdanford
 */
public class RegionParser implements Mapper<String,Region> {
    
    private static Pattern regPatt;  //Region pattern
    private static Pattern sregPatt;  //StrandedRegion pattern
    
    static { 
        regPatt = Pattern.compile("(\\w+):(\\d+)-(\\d+)");
        sregPatt = Pattern.compile("(\\w+):(\\d+)-(\\d+):(\\w+)");
    }
    
    private Genome genome;
    private int chromIndex, startIndex, endIndex, nameIndex, minLength;

    public RegionParser(Genome g) {
        genome = g;
        chromIndex = 0;
        startIndex = 1;
        endIndex = 2;
        nameIndex = 3;
        minLength = (Math.max(chromIndex, Math.max(startIndex, endIndex))) + 1;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(null)
     */
    public Region execute(String input) {
        String[] array = input.split("\\s+");
        String chrom = array[chromIndex];
        chrom = chrom.replaceFirst("chr", "");
        
        Matcher m = regPatt.matcher(chrom);
        if(m.matches()) { 
            chrom = m.group(1);
            int start = Integer.parseInt(m.group(2));
            int end = Integer.parseInt(m.group(3));
            return new Region(genome, chrom, start, end);
        } else { 
        	Matcher sm = sregPatt.matcher(chrom);
        	if(sm.matches()) { 
                chrom = sm.group(1);
                int start = Integer.parseInt(sm.group(2));
                int end = Integer.parseInt(sm.group(3));
                return new Region(genome, chrom, start, end);
        	}else{
	            if(array.length >= minLength) { 
	                int start = Integer.parseInt(array[startIndex]);
	                int end = Integer.parseInt(array[endIndex]);
	                if(nameIndex < array.length) {
	                    return new NamedRegion(genome, chrom, start, end, array[nameIndex]);
	                } else { 
	                    return new Region(genome, chrom, start, end);
	                }
	            } else { 
	                System.err.println("Line \"" + input + "\" doesn't have the correct length (" + minLength + ")");
	                return null;
	            }
        	}
        }
    }

}
