package edu.psu.compbio.seqcode.gse.datasets.alignments;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class Indel extends Region {
    public static final int INSERTION = 0, DELETION = 1, INVERSION = 2, REARRANGEMENT = 3;
    public static final String[] types = {"Insertion","Deletion","Inversion","Rearrangement"};

    private int type, size;

    public Indel(Genome g,
                 String chrom,
                 int start,
                 int end,
                 int size,
                 int type) {
        super(g,chrom,start,end);
        this.type = type;
        this.size = size;
    }
    public int getSize() {return size;}
    public int getType() {return type;}
    public String getTypeName() {return types[type];}
    public String toString() {
        return super.toString() + ", size=" + size + ", type=" + types[type];
    }

}