package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.*;

public class Seq {
	public int binsize;
	public boolean smoothing;
	public String orientation;
	public int length;
	public String chr;
	public int startcoord;
	public String sequence;
	
	public Seq(int startcoord, String chr, int binsize, boolean smoothing, String orientation, int length, String sequence ) {
		this.startcoord = startcoord;
		this.chr = chr;
		this.binsize = binsize;
		this.smoothing = smoothing;
		this.orientation = orientation;
		this.length = length;
		this.sequence = sequence;		
	}
}
