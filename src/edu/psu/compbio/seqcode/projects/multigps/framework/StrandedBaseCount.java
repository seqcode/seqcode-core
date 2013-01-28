package edu.psu.compbio.seqcode.projects.multigps.framework;

/**
 * StrandedBaseCount represents the sum of hit weights at a stranded base position in the genome.
 * The same class was called StrandedBase previously (renaming to reduce confusion)
 * Coordinate is the 5' end of the read. 
 * We do not store chromomosome info here because this class is always used
 * 		in the context of a chromosome or a region, no need to differentiate.
 * It records the number of reads mapped to this base position. 
 * For deeply-seq dataset, the count is typically higher.
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
		
}
