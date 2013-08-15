package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.*;

public class Vec {
	public String orientation;
	public int range;
	public String chr;
	public Map<Integer,Integer> tags = new TreeMap<Integer, Integer>();
	public int binsize;
	public boolean smoothing;
	
	public Vec(int range, String chr, String orientation, boolean smoothing, int binsize, Map<Integer, Integer> tags) {
		this.range = range;
		this.chr = chr;
		this.orientation = orientation;
		this.smoothing = smoothing;
		this.binsize = binsize;
		this.tags = tags;
	}
	 public Vec getSub(int midpoint, int range){
		 int start = midpoint-range;
		 int end = midpoint+range;
		 Map<Integer, Integer> rettags = new TreeMap<Integer, Integer>();
		 for(int ke: this.tags.keySet()){
			 if(ke>=start && ke < end){
				 rettags.put(ke, this.tags.get(ke));
			 }
		 }
		 Vec ret = new Vec(this.range,this.chr,this.orientation,this.smoothing,this.binsize,rettags);
		 return ret;
	}
}
