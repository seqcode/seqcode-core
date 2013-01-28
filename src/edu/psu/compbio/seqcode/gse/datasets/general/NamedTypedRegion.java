package edu.psu.compbio.seqcode.gse.datasets.general;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class NamedTypedRegion extends NamedStrandedRegion implements Typed, Stranded {

    private char strand;
    private String type;

    public NamedTypedRegion (Genome g, String c, int start, int end, String name, String type, char strand) {
        super(g,c,start,end,name,strand);
        this.strand = strand;
        this.type = type;
    }

    public String getType() {return type;}
    public char getStrand() {return strand;}
}
