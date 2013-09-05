package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;

public class ChexmixSandbox {
	
	public static int[] getCompositeFromBlLisr(Map<BindingLocation, String> blmap){
		List<BindingLocation> listOfBl = new ArrayList<BindingLocation>(blmap.keySet());
		int[] ret  = new int[listOfBl.get(0).vecpos.tags.size()*2];
		
		for(int i=0; i<listOfBl.size(); i++){
			List<Integer> keylist = new ArrayList<Integer>(listOfBl.get(i).vecpos.tags.keySet());
			if(blmap.get(listOfBl.get(i)).equals("+")){
				for(int j=0; j<keylist.size(); j++){
					if(i==0){
						ret[j] = listOfBl.get(i).vecpos.tags.get(keylist.get(j));
					}
					else{
						ret[j] = ret[j] + listOfBl.get(i).vecpos.tags.get(keylist.get(j));
					}
					
				}
				for(int j=keylist.size()-1; j>=0; j--){
					if(i==0){
						ret[2*keylist.size()-1-j] = listOfBl.get(i).vecneg.tags.get(keylist.get(j));
					}
					else{
						ret[2*keylist.size()-1-j] = ret[2*keylist.size()-1-j] + listOfBl.get(i).vecneg.tags.get(keylist.get(j));
					}
				}
			}
			else{
				for(int j=0; j<keylist.size(); j++){
					if(i==0){
						ret[keylist.size()+j] = listOfBl.get(i).vecpos.tags.get(keylist.get(j));
					}
					else{
						ret[keylist.size()+j] = ret[keylist.size()+j] + listOfBl.get(i).vecpos.tags.get(keylist.get(j));
					}
					
				}
				for(int j=keylist.size()-1; j>=0; j--){
					if(i==0){
						ret[keylist.size()-1-j] = listOfBl.get(i).vecneg.tags.get(keylist.get(j));
					}
					else{
						ret[keylist.size()-1-j] = ret[keylist.size()-1-j] + listOfBl.get(i).vecneg.tags.get(keylist.get(j));
					}
				}
				
			}
			
			
			
		}
		return ret;
	}

}
