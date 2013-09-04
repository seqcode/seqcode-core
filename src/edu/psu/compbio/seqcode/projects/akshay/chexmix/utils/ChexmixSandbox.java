package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;

public class ChexmixSandbox {
	
	public static int[] getCompositeFromBlLisr(List<BindingLocation> bllist){
		int[] ret  = new int[bllist.get(0).vecpos.tags.size()*2];
		
		for(int i=0; i<bllist.size(); i++){
			List<Integer> keylist = new ArrayList<Integer>(bllist.get(i).vecpos.tags.keySet());
			for(int j=0; j<keylist.size(); j++){
				if(i==0){
					ret[j] = bllist.get(i).vecpos.tags.get(keylist.get(j));
				}
				else{
					ret[j] = ret[j] + bllist.get(i).vecpos.tags.get(keylist.get(j));
				}
				
			}
			for(int j=keylist.size()-1; j>=0; j--){
				if(i==0){
					ret[2*keylist.size()-1-j] = bllist.get(i).vecneg.tags.get(keylist.get(j));
				}
				else{
					ret[2*keylist.size()-1-j] = ret[2*keylist.size()-1-j] + bllist.get(i).vecneg.tags.get(keylist.get(j));
				}
			}
			
		}
		return ret;
	}

}
