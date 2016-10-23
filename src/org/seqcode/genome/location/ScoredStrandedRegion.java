package org.seqcode.genome.location;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import org.seqcode.genome.Genome;


public class ScoredStrandedRegion extends ScoredRegion implements Stranded {
    private char strand;
    
    public ScoredStrandedRegion(ScoredStrandedRegion copied) {
        super(copied);
        strand = copied.strand;
    }

    public ScoredStrandedRegion(Genome g, String c, int start, int end, double score, char strand) {
        super(g,c,start,end,score);
        this.strand = strand;
    }
    
    public ScoredStrandedRegion(Region r, double score, char strand) {
        super(r.getGenome(), r.getChrom(),r.getStart(), r.getEnd(),score);
        this.strand = strand;
    }
    
    public ScoredStrandedRegion(Genome g, DataInputStream dis) throws IOException { 
        super(g,dis);
        strand = dis.readChar();
    }
    
    public void save(DataOutputStream dos) throws IOException { 
        super.save(dos);
        dos.writeChar(strand);
    }

    public char getStrand() {return strand;}

    public boolean equals(Object o) {
        if (o instanceof ScoredStrandedRegion) {
            ScoredStrandedRegion other = (ScoredStrandedRegion)o;
            return super.equals(other) && other.strand == strand;
        } else {
            return false;
        }
    }
    
    public ScoredStrandedRegion expand(int upstream, int downstream) {
        if (strand == '+') {
            int ns = getStart() - upstream;
            int ne = getEnd() + downstream;
            if (ns < 1) {ns = 1;}
            return new ScoredStrandedRegion(getGenome(),getChrom(),ns,ne,score, strand);
        } else if (strand == '-') {
            int ns = getStart() - downstream;
            int ne = getEnd() + upstream;
            if (ns < 1) {ns = 1;}
            return new ScoredStrandedRegion(getGenome(),getChrom(),ns,ne,score, strand);                
        } else {
            throw new IllegalArgumentException("Strand isn't + or - so I don't know what to do");
        }

    }
    public String toString() { 
        String str = getLocationString() + ":" + strand;
        str += " (" + getScore() + ")";
        return str;
    }
    public String toTabString() { 
        String str = getLocationString() + ":" + strand;
        str += "\t" + getScore();
        return str;
    }
}
