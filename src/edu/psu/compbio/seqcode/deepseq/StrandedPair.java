package edu.psu.compbio.seqcode.deepseq;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;

/**
 * StrandedPair represents a pair of reads that are paired together.
 * Coordinate is the 5' end of each hit. 
 * Weight is determined by implementation of hit loading. 
 * 
 * @author mahony
 *
 */
public class StrandedPair implements Comparable<StrandedPair>{
	private Genome gen;
	private int r1Chrom, r2Chrom; // int codes for chroms - convert with Genome
	private char r1Strand, r2Strand; 
	private int r1Coordinate, r2Coordinate;
	private boolean sameChrom;
	private float weight;
	
	public StrandedPair(Genome g, int r1Chr, int r1Coord, char r1Str, int r2Chr, int r2Coord, char r2Str, boolean sameChr, float w){
		gen=g;
		r1Chrom = r1Chr;
		r2Chrom = r2Chr;
		r1Strand = r1Str;
		r2Strand = r2Str;
		r1Coordinate = r1Coord;
		r2Coordinate = r2Coord;
		sameChrom = sameChr;
		weight = w;
	}
	
	/**
	 * Get the pair midpoint. Returns null if reads on different strands.
	 * @return
	 */
	public Point getMidpoint(){
		if(!sameChrom)
			return null;
		else{
			return new Point(gen, gen.getChromName(r1Chrom), (r1Coordinate+r2Coordinate)/2);
		}
	}
	
	/** 
	 * Returns distance between read 5' positions if this is a proper concurrent pair. 
	 * Returns -1 otherwise. 
	 * @return
	 */
	public int getFragmentSize(){
		if(!sameChrom)
			return -1;
		else if((r1Coordinate<r2Coordinate && r1Strand=='+' && r2Strand=='-') || (r2Coordinate<r1Coordinate && r2Strand=='+' && r1Strand=='-')){
			return Math.abs(r2Coordinate-r1Coordinate);
		}else{
			return -1;
		}
	}

	public void setR1Strand(char strand) {
		this.r1Strand = strand;
	}
	public void setR2Strand(char strand) {
		this.r2Strand = strand;
	}

	public char getR1Strand() {
		return r1Strand;
	}
	public char getR2Strand() {
		return r2Strand;
	}

	public void setR1Coordinate(int coordinate) {
		this.r1Coordinate = coordinate;
	}
	public void setR2Coordinate(int coordinate) {
		this.r2Coordinate = coordinate;
	}

	public int getR1Coordinate() {
		return r1Coordinate;
	}
	public int getR2Coordinate() {
		return r2Coordinate;
	}
	
	public void setR1Chrom(int chrom){
		if(r1Chrom == r2Chrom)
			sameChrom=true;
		this.r1Chrom=chrom;
	}
	public void setR2Chrom(int chrom){
		if(r1Chrom == r2Chrom)
			sameChrom=true;
		this.r2Chrom=chrom;
	}
	public String getR1Chrom(){
		return gen.getChromName(r1Chrom);
	}
	public String getR2Chrom(){
		return gen.getChromName(r2Chrom);
	}
	
	public boolean pairFromSameChrom(){ return sameChrom; }

	public void setWeight(float weight) {
		this.weight = weight;
	}
	public float getWeight() {
		return weight;
	}
	
	// sort according to R1 coordinate, then by R2 coordinate, then by strands
	public int compareTo(StrandedPair b) {
		int result;
        if (r1Chrom==b.r1Chrom) {   //R1 coord
            result = r1Coordinate - b.r1Coordinate;
        } else {
            result = r1Chrom > b.r1Chrom ? +1 : r1Chrom < b.r1Chrom ? -1 : 0;
        }
        if (result == 0) {   //R2 coord
            if (r2Chrom == b.r2Chrom) {
                result = r2Coordinate - b.r2Coordinate;
            } else {
            	result = r2Chrom > b.r2Chrom ? +1 : r2Chrom < b.r2Chrom ? -1 : 0;
            }
        }
        if (result == 0)    //R1 strand
            result = (r1Strand==b.r1Strand ? 0 : r1Strand=='+'?1:-1);
        if (result == 0) {    //R2 strand
        	result = (r2Strand==b.r2Strand ? 0 : r2Strand=='+'?1:-1);
        }
        return result;
	}
			
}
