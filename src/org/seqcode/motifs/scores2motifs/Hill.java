package org.seqcode.motifs.scores2motifs;

import org.seqcode.genome.location.Region;

public class Hill {

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
}
