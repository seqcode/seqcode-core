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
	
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		
		Seq se = (Seq) obj;
		return this.orientation == se.orientation && this.range == se.range && this.chr == se.chr && this.midpoint == se.midpoint;
	}
	
	@Override
	public int hashCode(){
		int result = 17;
		int code = (int) this.range;
		code+= (int) this.midpoint;
		code += (int) (this.chr == null ? 0 :this.chr.hashCode());
		code+= (int) (this.orientation == null ? 0 : this.orientation.hashCode());
		
		result = result*37 + code;
		return result;
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
	
	public static void main(String[] args){
		
		Map<Integer,Integer> test = new HashMap<Integer, Integer>();
		System.out.println(test.size());
		test.put(3, 4);
		System.out.println(test.size());
		test.remove(3);
		System.out.println(test.size());
	}
}
