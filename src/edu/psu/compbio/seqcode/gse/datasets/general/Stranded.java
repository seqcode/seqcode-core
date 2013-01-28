package edu.psu.compbio.seqcode.gse.datasets.general;

/**
 * Interface for anything (generally a Region) that has strand information.
 */
public interface Stranded {
    public char getStrand();
}
