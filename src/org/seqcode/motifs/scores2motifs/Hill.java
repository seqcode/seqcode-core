package org.seqcode.motifs.scores2motifs;

import org.seqcode.genome.location.Region;

public class Hill implements Comparable<Hill>{

	protected Region hillLoc;
	protected Double score;
	protected String seq;
	protected int[] kmerProfile=null;
	
	public Hill(Region r, Double s, String seq){
		hillLoc = r;
		score = s;
		this.seq=seq;
	}
	
	//Settor
	public void setKmerProfile(int[] kp){ kmerProfile=kp; }
	
	//Gettors
	public Region getHillRegion(){return hillLoc;}
	public Double getScore(){return score;}
	public String getSeq(){return seq;}
	public int[] getKmerProfile(){return kmerProfile;}

	
	public int compareTo(Hill o) {
		return o.score.compareTo(this.score);
	}
}
