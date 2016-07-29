package org.seqcode.genome.sequence;

/**
 * Sequence (string) with an attached score
 * @author mahony
 *
 */
public class ScoredSequence {
	private String seq;
	private double score;
	
	public ScoredSequence(String seq, double score){
		this.seq = seq;
		this.score = score;
	}
	
	//Accessors
	public String getSeq(){return seq;}
	public double getScore(){return score;}
	public void setSeq(String s){seq=s;}
	public void setScore(double s){score=s;}
}