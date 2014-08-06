package edu.psu.compbio.seqcode.genome.location;

/**
 * Interface for anything (generally a Region) that has strand information.
 */
public interface Stranded {
    public char getStrand();
}
