package org.seqcode.deepseq;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.ScoredStrandedPoint;
import org.seqcode.genome.location.ScoredStrandedRegion;

/**
 * StrandedBaseCount represents the sum of hit weights at a stranded base position in the genome.
 * The same class was called StrandedBase previously (renaming to reduce confusion)
 * Coordinate is the 5' end of the read. 
 * We do not store chromomosome info here because this class is always used
 * 		in the context of a chromosome or a region, no need to differentiate.
 * It records the number of reads mapped to this base position. 
 * 
 * @author yguo
 *
 */
public class StrandedBaseCount implements Comparable<StrandedBaseCount>{
	private char strand; 
	private int coordinate;
	private float count;
	
	public StrandedBaseCount(char strand, int coord, float count){
		this.setStrand(strand);
		this.setCoordinate(coord);
		this.setCount(count);
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}

	public char getStrand() {
		return strand;
	}

	public void setCoordinate(int coordinate) {
		this.coordinate = coordinate;
	}

	public int getCoordinate() {
		return coordinate;
	}

	public void setCount(float count) {
		this.count = count;
	}

	public float getCount() {
		return count;
	}
	// sort according to coordinate, considering strand
	public int compareTo(StrandedBaseCount b) {
		double diff = coordinate-b.coordinate;
		return diff==0?
				(strand==b.strand?0:strand=='+'?1:-1):	// same coord, compare strand
					(diff<0?-1:1);						// diff coord
	}
	
	public String toString(){
		return coordinate+" "+strand+" "+count;
	}
		
	//Makes this into a regular Region-derived object
	public ScoredStrandedPoint toScoredStrandedPoint(Genome g, String chr){
		return new ScoredStrandedPoint(g, chr, coordinate, count, strand);
	}
	
	//Makes the pseudo extended read into a regular Region-derived object
	public ScoredStrandedRegion expandToScoredStrandedRegion(Genome g, String chr, int upstreamext, int downstreamext){
		return (new ScoredStrandedPoint(g, chr, coordinate, count, strand).expand(upstreamext, downstreamext));
	}
}
