package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import edu.psu.compbio.seqcode.genome.location.Region;

public class AlignBlock{
	protected int length;
	protected int refStart;
	
	public AlignBlock(int len, int refStart){
		length = len;
		this.refStart = refStart;		
	}
	
	//Accessors
	public int getLength(){return length;}
	public int getReferenceStart(){return refStart;}
	public int getReferenceEnd(){return refStart+length-1;} //Assumes contiguous alignment
	
	/**
	 * Checks for overlap, assuming that they are on the same chromosome
	 */
	public boolean overlaps(Region r) {
	    int refEnd = refStart+length-1;
		if (refStart <= r.getStart() && refEnd >= r.getStart()) {
	      return true;
	    }
	    if (r.getStart() <= refStart && r.getEnd() >= refStart) {
	      return true;
	    }
	    return false;
	  }
}
