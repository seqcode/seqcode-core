package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis.LoadTags;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;

public class QueryTags {
	public int midpoint;
	public String chr;
	public int range;
	public int Smoothsize;
	public Map<Integer, Integer> tags = new TreeMap<Integer, Integer>();
	
	public QueryTags(int midpoint, int range, String chr, int Smoothsize) {
		this.midpoint = midpoint;
		this.range = range;
		this.chr = chr;
	}
	
	public Vec getTags(LoadTags loader, String orientation){
		Vec ret = null;
		Map<Integer,Integer> tags = new TreeMap<Integer,Integer>();
		String[] tmp = this.chr.split("\\.");
		String ch=tmp[0].replaceFirst("chr", "");
		ch=ch.replaceFirst("^>", "");
		int chrID = loader.chrom2ID.get(ch);
		int j = (orientation == "+")?0:1;
		if(loader.fivePrimePos[chrID][j]!= null){
			int[] tempStarts = loader.fivePrimePos[chrID][j];
			if(tempStarts.length != 0){
				int start_ind = Arrays.binarySearch(tempStarts, this.midpoint-this.range);
				//debug line
				System.out.println(start_ind);
				System.out.println(tempStarts.length);
				System.out.println(tempStarts[tempStarts.length-1]);
				//end
				if( start_ind < 0 ) { start_ind = -start_ind - 1; }
				for(int k=this.midpoint-this.range; k<this.midpoint+this.range; k++){
					if(k == tempStarts[start_ind]){
						tags.put(k, (int)loader.fivePrimeCounts[chrID][j][start_ind]);
						start_ind++;
					}
					else{
						tags.put(k, 0);
					}
				}
				
			}
		}
		ret = new Vec(this.range, this.midpoint, this.chr, orientation, this.Smoothsize, 0,tags);
		return ret;
	}
	
	public static void main(String[] args){
		String strand="+";
		int pp = (strand == "-"?0:1);
		System.out.println(pp);
	}
}
		

