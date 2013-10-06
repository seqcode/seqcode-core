package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.CustomReturn;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Membership;

public class ChexmixSandbox {
	
	public static int[] getCompositeFromBlLisr(List<BindingLocation>listofbl, Map<BindingLocation, String> blmap){
		
		int[] ret  = new int[listofbl.get(0).vecpos.tags.size()*2];
		
		for(int i=0; i<listofbl.size(); i++){
			List<Integer> keylist = new ArrayList<Integer>(listofbl.get(i).vecpos.tags.keySet());
			if(blmap.get(listofbl.get(i)).equals("+")){
				for(int j=0; j<keylist.size(); j++){
					if(i==0){
						ret[j] = listofbl.get(i).vecpos.tags.get(keylist.get(j));
					}
					else{
						ret[j] = ret[j] + listofbl.get(i).vecpos.tags.get(keylist.get(j));
					}
					
				}
				for(int j=keylist.size()-1; j>=0; j--){
					if(i==0){
						ret[2*keylist.size()-1-j] = listofbl.get(i).vecneg.tags.get(keylist.get(j));
					}
					else{
						ret[2*keylist.size()-1-j] = ret[2*keylist.size()-1-j] + listofbl.get(i).vecneg.tags.get(keylist.get(j));
					}
				}
			}
			else{
				for(int j=0; j<keylist.size(); j++){
					if(i==0){
						ret[keylist.size()+j] = listofbl.get(i).vecpos.tags.get(keylist.get(j));
					}
					else{
						ret[keylist.size()+j] = ret[keylist.size()+j] + listofbl.get(i).vecpos.tags.get(keylist.get(j));
					}
					
				}
				for(int j=keylist.size()-1; j>=0; j--){
					if(i==0){
						ret[keylist.size()-1-j] = listofbl.get(i).vecneg.tags.get(keylist.get(j));
					}
					else{
						ret[keylist.size()-1-j] = ret[keylist.size()-1-j] + listofbl.get(i).vecneg.tags.get(keylist.get(j));
					}
				}
				
			}
			
			
			
		}
		return ret;
	}
	
	
	public static List<Membership> perfromReAssignmentOfMembership(List<Membership> initial_membership, List<int[]> seed_profiles){
		List<Membership> ret = initial_membership;
		for(int i=0; i< ret.size(); i++){
			double max_pcc = ret.get(i).getPCC();
			int k=0;
			while(k< seed_profiles.size()){
				if(k+1 != ret.get(i).membership){
					CustomReturn temp_cr = ret.get(i).bl.scanConcVecWithBl(seed_profiles.get(k), ret.get(i).cr.maxvec.range);
					if(temp_cr.pcc > max_pcc){
						max_pcc = temp_cr.pcc;
						ret.get(i).cr = temp_cr;
						ret.get(i).membership = k+1;
					}
				}
			}
		}
		return ret;
	}
	
	public static double[][] getClusterAssesmsent(List<Membership> assignment, List<int[]> seed_profiles){
		// return is a 2-d array.. with rows as the seed profile and column as the clusters
		double[][] ret = new double[seed_profiles.size()][seed_profiles.size()];
		for(int i=0; i<seed_profiles.size(); i++){
			for(int j=0; j < seed_profiles.size(); j++){
				if(i==j){
					double add = 0.0;
					int count = 0;
					for(int k=0; k < assignment.size(); k++){
						if(assignment.get(k).membership == j+1){
							add = add + assignment.get(k).getPCC();
							count++;
						}
					}
					ret[i][j] = add/count;
				}
				else{
					double add = 0.0;
					int count = 0;
					for(int k=0; k< assignment.size(); k++){
						if(assignment.get(k).membership == j+1){
							add = add + assignment.get(k).bl.scanConcVecWithBl(seed_profiles.get(i), assignment.get(k).cr.maxvec.range).pcc;
							count++;
						}
					}
					ret[i][j] = add/count;
				}
			}
		}
		return ret;
	}

}
