package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.utils.*;

public class Vec {
	public String orientation;
	public int range;
	public String chr;
	public Map<Integer,Integer> tags = new TreeMap<Integer, Integer>();
	public int binsize;
	public boolean smoothing;
	public int midpoint;
	
	public Vec(int range, int midpoint, String chr, String orientation, boolean smoothing, int binsize, Map<Integer, Integer> tags) {
		this.range = range;
		this.chr = chr;
		this.orientation = orientation;
		this.smoothing = smoothing;
		this.binsize = binsize;
		this.tags = tags;
		this.midpoint = midpoint;
	}
	
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		
		Vec ve = (Vec) obj;
		return this.orientation == ve.orientation && this.range == ve.range && this.chr == ve.chr && this.binsize == ve.binsize && this.smoothing == ve.smoothing;
	}
	
	@Override
	public int hashCode(){
		int result = 17;
		int code = (int) this.range;
		code+= (int) this.midpoint;
		code += (int) (this.chr == null ? 0 :this.chr.hashCode());
		code+= (this.smoothing ? 1231 : 1237);
		code+= (int) (this.orientation == null ? 0 : this.orientation.hashCode());
		
		result = result*37 + code;
		return result;
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
		 Vec ret = new Vec(range, midpoint, this.chr,this.orientation,this.smoothing,this.binsize,rettags);
		 return ret;
	}
	 
}
