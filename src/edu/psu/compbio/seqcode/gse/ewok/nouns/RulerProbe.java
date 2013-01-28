package edu.psu.compbio.seqcode.gse.ewok.nouns;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class RulerProbe extends Region {

    private double intensity;
    private int probeid;
    private char probeSide;

    public RulerProbe (Genome g,
                       String c,
                       int start, 
                       int end,
                       int probeid,
                       char probeSide,
                       double intensity) {
        super(g,c,start,end);
        this.intensity = intensity;
        this.probeid = probeid;
        this.probeSide = probeSide;
    }

    public double getIntensity() {return intensity;}
    public int getProbeID() {return probeid;}
    public char getProbeSide() {return probeSide;}                       
}
