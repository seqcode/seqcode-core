package org.seqcode.projects.akshay.chexmix.utils;

import java.io.File;
import java.util.*;

import org.seqcode.projects.akshay.chexmix.datasets.*;


public class Smoothing {
	public Vec doSmoothing(Vec inputvec, int smoothsize){
		Map<Integer, Integer> smoothtags = new TreeMap<Integer, Integer>();
		if(smoothsize == 3){
			List<Integer> keylist = new ArrayList<Integer>(inputvec.tags.keySet());
			smoothtags.put(keylist.get(0), inputvec.tags.get(keylist.get(0)));
			smoothtags.put(keylist.get(keylist.size()-1), inputvec.tags.get(keylist.get(keylist.size()-1)));
			for(int i=1; i<keylist.size()-1; i++){
				smoothtags.put(keylist.get(i),(inputvec.tags.get(keylist.get(i))+inputvec.tags.get(keylist.get(i-1))+inputvec.tags.get(keylist.get(i+1)))/3);
			}
		}
		if(smoothsize == 5){
			List<Integer> keylist = new ArrayList<Integer>(inputvec.tags.keySet());
			smoothtags.put(keylist.get(0), inputvec.tags.get(keylist.get(0)));
			smoothtags.put(keylist.get(1), inputvec.tags.get(keylist.get(1)));
			smoothtags.put(keylist.get(keylist.size()-1), inputvec.tags.get(keylist.get(keylist.size()-1)));
			smoothtags.put(keylist.get(keylist.size()-2), inputvec.tags.get(keylist.get(keylist.size()-2)));
			for(int i=2; i<keylist.size()-2; i++){
				smoothtags.put(keylist.get(i),(3*inputvec.tags.get(keylist.get(i))+2*inputvec.tags.get(keylist.get(i-1))+2*inputvec.tags.get(keylist.get(i+1))+inputvec.tags.get(keylist.get(i-2))+inputvec.tags.get(keylist.get(i+2)))/9);
			}
		}
		Vec ret = new Vec(inputvec.range, inputvec.midpoint, inputvec.chr, inputvec.orientation, smoothsize, inputvec.binsize, smoothtags);
		return ret;
		
	}
	
	public static void main(String[] args){
		String chr = "chr7";
		String[] tmp = chr.split("\\.");
		System.out.println(tmp[0]);
		chr=tmp[0].replaceFirst("chr", "");
		System.out.println(chr);
		chr=chr.replaceFirst("^>", "");
		System.out.println(chr);
		int a=3;
		int b=4;
		int c=5;
		int d= (1)/5;
		System.out.println(d); 
		
	}
}
