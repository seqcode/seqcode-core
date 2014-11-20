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
public class HitPair{
	public int r1Pos;
	public String r2Chr;
	public int r2Pos;
	public int r2Strand; // 0 for '-', 1 for '+'
	public float pairWeight;
	public HitPair(int p1, String c2, int p2, int s2, float w){
		r1Pos = p1; r2Chr = c2; r2Pos = p2; r2Strand = s2; pairWeight=w;
	}
}
