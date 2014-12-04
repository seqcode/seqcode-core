package edu.psu.compbio.seqcode.deepseq;

/**
 * Mini class to hold four bits of information on the hit pair:
 * (assumes that R1 hit chromosome and strand are stored in the data structure itself
 *  - R1 hit 5' pos
 *  - R2 hit chr
 *  - R2 hit 5' pos
 *  - R2 hit strand
 *   
 * @author mahony
 */
public class HitPair implements Comparable<HitPair>{
	public int r1Pos;
	public String r2Chr;
	public int r2Pos;
	public int r2Strand; // 0 for '-', 1 for '+'
	public float pairWeight;
	public HitPair(int p1, String c2, int p2, int s2, float w){
		r1Pos = p1; r2Chr = c2; r2Pos = p2; r2Strand = s2; pairWeight=w;
	}
	
	// sort according to R1 coordinate, then by R2 coordinate, then by strands
	public int compareTo(HitPair b) {
		int result = r1Pos - b.r1Pos;
        
		if (result == 0) {   //R2 coord
            if (r2Chr.equals(b.r2Chr)) {
                result = r2Pos - b.r2Pos;
            } else {
            	result = r2Chr.compareTo(b.r2Chr);
            }
        }
        if (result == 0) {    //R2 strand
        	result = (r2Strand==b.r2Strand ? 0 : r2Strand=='+'?1:-1);
        }
        return result;
	}

}
