package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.*;

public class Seq {
	public String orientation;
	public int range;
	public String chr;
	public String sequence;
	public int midpoint;
	
	public Seq(int midpoint, int range, String chr, String orientation, String sequence ) {
		this.range = range;
		this.chr = chr;
		this.orientation = orientation;
		this.sequence = sequence;
		this.midpoint = midpoint;
	}
	
	public Seq getSub(int midpoint, int range){
		Seq ret=null;
		if(midpoint-range/2 < this.midpoint-this.range/2 && midpoint+range/2 >this.midpoint + this.range/2 ){
			System.err.println("getSub err: the sequence requested is not a subsequence of the original sequence!");
		}
		else{
			ret = new Seq(midpoint, range, this.chr, this.orientation, this.sequence.substring(midpoint-range/2, midpoint+range/2));
		}
		return ret;
	}
}
