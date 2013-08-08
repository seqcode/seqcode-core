package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.*;

public class BindingLocation {

	private Seq seqpos;
	private Seq seqneg;
	private Vec vecpos;
	private Vec vecneg;
	private String chr;
	private List<Integer> coords = new ArrayList<Integer>();
	private Map<Integer,Integer> tags = new TreeMap<Integer,Integer>();
	
	public BindingLocation(int midpoint, String chr, int length) {
		this.chr = chr;
		int start = midpoint - length/2;
		int end = midpoint + length/2;
		this.coords.add(start);
		this.coords.add(end);
	}
	
	public void fillSeqs(){
		
	}
	
	
}
