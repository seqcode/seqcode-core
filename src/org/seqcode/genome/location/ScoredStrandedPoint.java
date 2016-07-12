package org.seqcode.genome.location;

import org.seqcode.genome.Genome;

public class ScoredStrandedPoint extends ScoredPoint implements Stranded {
	
	private char strand;
	
    public ScoredStrandedPoint (Genome g, String c, int position, double s, char str) {
    	super(g, c, position, s);
    	strand = str;
    }
    
    public ScoredStrandedPoint (Point p, double s, char str) {
    	super(p.getGenome(), p.getChrom(), p.getLocation(), s);
    	strand = str;
    }
    
    public char getStrand() { return strand; }
    
    public ScoredStrandedPoint clone() {
        return new ScoredStrandedPoint(getGenome(), getChrom(), getLocation(), getScore(), strand);
    }
    
    public String toString() {
        return String.format("%s:%c (%.3f)", getLocationString(), strand, getScore());
    }
    
    public boolean equals(Object o) {
        if (o instanceof ScoredStrandedPoint) {
            ScoredStrandedPoint r = (ScoredStrandedPoint)o;
            if(!super.equals(r)) { return false; }
            if(strand != r.strand) { return false; }
            return true;
        } else {
            return false;
        }
    }
}
