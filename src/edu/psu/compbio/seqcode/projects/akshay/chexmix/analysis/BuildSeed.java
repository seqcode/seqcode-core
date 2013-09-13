package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.utils.Pearson;


public class BuildSeed {
	
	/**
	 * ALLPAIRS is a hashmap with keys as the BLpair object (edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BLpair)
	 * the value of this hashmap is a CustomReturn object. The value is obtained by scanning the two BL's in BLpair using the scanWithBl method in BindingLocation class
	 * filled by the constructor
	 */
	
	public Map<BLpair,CustomReturn> allpairs = new HashMap<BLpair, CustomReturn>();
	/**
	 * REMOVED is hashmap that contains all the binding locations that are far away from the center (based on the cutoff). 
	 * The value is just a dummy integer (1 in this case)
	 * Filled be the updateRemoved method
	 */
	public Map<BindingLocation,Integer> removed = new HashMap<BindingLocation,Integer>();
	/**
	 * PCCPAIRWISE is a hashmap whose keys are all the binding locations given and whose value is the sum of PCC from every other binding location in the list.
	 * Also filled by the getCenter method
	 */
	public Map<BindingLocation,Double> pccpairwise = new HashMap<BindingLocation,Double>();
	/**
	 * SEED is an ArrayList that contains all the bindinglocations that have been chosen to construct the seed profile. Filled by the getSeedList method
	 */
	public List<BindingLocation> seed;
	/**
	 * CENTER is the bindinglocation that is closest to every other binding locations. This is the filled by the getCenter method
	 */
	public BindingLocation center;
	
	
	/**
	 * Initiats the class object that should be used to construct the seed profile
	 * @param topBlList (The entire list of bindinglocations that are used to construct the seed profile)
	 * @param range (The size of the range; this parameter is taken from the config object)
	 * @param smoothsize (Smoothing bin size; this parameter is taken from the config object)
	 */
	public BuildSeed(List<BindingLocation> topBlList, Config conf) {
		for(int i=0; i<topBlList.size()-1; i++){
			for(int j=i+1; j< topBlList.size(); j++ ){
				if(i != j){
					//debugging line
					//System.out.println(i+"to"+j);
					BLpair pair = new BLpair(topBlList.get(i),topBlList.get(j));
					CustomReturn pccmax = topBlList.get(i).scanBlWithBl(topBlList.get(j), conf.getIntSize());
					//debug line
					//System.out.println(pccmax.pcc);
					allpairs.put(pair, pccmax);
				}
			}
		}
		if(conf.useCenter()){
			this.getCenter();
			List<BindingLocation> ret = new ArrayList<BindingLocation>();
			this.updateRemoved(center,conf.getSeedCutoff());
			for(BindingLocation bl : pccpairwise.keySet()){
				if(!removed.containsKey(bl)){
					ret.add(bl);
				}
			}
			this.seed = ret;
		}
		else{
			this.seed = new ArrayList<BindingLocation>(pccpairwise.keySet());
		}
		
	}
	
	/**
	 * Used by the constructor; should not be called otherwise
	 * Gets the bindinglocation that is closest to every other bindinglocation (By maximising the pair-wise pearson correlation co-eff)
	 * @return
	 */
	private void getCenter(){
		for (BLpair pair: allpairs.keySet()){
			if(pccpairwise.containsKey(pair.BL1)){
				pccpairwise.put(pair.BL1, pccpairwise.get(pair.BL1) + allpairs.get(pair).pcc);
			}
			else{
				pccpairwise.put(pair.BL1, allpairs.get(pair).pcc);
			}
			if(pccpairwise.containsKey(pair.BL2)){
				pccpairwise.put(pair.BL2, pccpairwise.get(pair.BL2)+allpairs.get(pair).pcc);
			}
			else{
				pccpairwise.put(pair.BL2, allpairs.get(pair).pcc);
			}
		}
		BindingLocation ret=null;
		double max=-1000;
		
		for(BindingLocation bl : pccpairwise.keySet()){
			if(pccpairwise.get(bl)>max){
				max =pccpairwise.get(bl);
				ret = bl;
			}
		}
		this.center = ret;
	}
	
	/**
	 * Used by the getSeedList; should not be used otherwise
	 * Given the center; this method marks all those binding locations that are far away from the center (based on the cutoff) by adding them 
	 * to the removed hashmap
	 * @param center
	 * @param cutoff
	 */
	private void updateRemoved(BindingLocation center, double cutoff){
		for(BLpair pair : allpairs.keySet()){
			if(pair.BL1 == center){
				if(allpairs.get(pair).pcc < cutoff){
					if(!removed.containsKey(pair.BL2)){
					removed.put(pair.BL2, 1);
					}
				}
			}
			if(pair.BL2 == center){
				if(allpairs.get(pair).pcc < cutoff){
					if(!removed.containsKey(pair.BL1)){
						removed.put(pair.BL1, 1);
					}
				}
			}
		}
		
	}
	
	/**
	 * Given a list of different seed profiles, this method finds the best one by doing a post assessment.
	 * It calculates the average PCC of the seed with every location that was used to build it. The profile with highest avg is returned.
	 * This function is generally called by one of the schemes and should not be called separately
	 * @param allseeds
	 * @return
	 */
	public int[] doPostAssessment(Map<BindingLocation, int[]> allseeds, Config conf ){
		double maxavgpcc=-1000.0;
		int[] bestalignment=null;
		for(BindingLocation bl: allseeds.keySet()){
			double avgpcc=0.0;
			int[] profile = allseeds.get(bl);
			for(BindingLocation jbl : this.seed){
				CustomReturn tempRet = jbl.scanConcVecWithBl(profile, conf.getIntSize());
				avgpcc = avgpcc+tempRet.pcc;
			}
			avgpcc = avgpcc/allseeds.size();
			if(avgpcc> maxavgpcc){
				maxavgpcc = avgpcc;
				bestalignment = profile;
			}
		}
		return bestalignment;
	}
	
	/**
	 * This is the scheme-1 for building the seed.
	 * This scheme works when the center (exemplar point) is the true representation of the data (other binding locations considered)
	 * This method generates n-1 profiles. Where n is the number in the seed list. Each profile is constructed based on a displacement in the center
	 *  that maximizes similarity between center and a particular bl. The best profile is selected by doing a post assessment.
	 * @return
	 */
	public int[] executeScheme1(Config conf){
		int[] bestprofile = null;
		Map<BindingLocation,int[]> allseeds = new HashMap<BindingLocation,int[]>();
		for(BindingLocation bl: seed){
			int[] profile=null;
			BLpair temppair = (allpairs.containsKey(new BLpair(center, bl)) ? new BLpair(center,bl): new BLpair(bl,center));
			CustomReturn tempPos = this.allpairs.get(temppair);
			List<Integer> tempaddToProfile = temppair.BL1.getConcatenatedTags(tempPos.maxVec1.midpoint, tempPos.maxVec1.range, tempPos.maxVec1.orientation); 
			profile = new int[tempaddToProfile.size()];
			for(int i=0; i<tempaddToProfile.size(); i++){
				profile[i] = tempaddToProfile.get(i);
			}
			tempaddToProfile = temppair.BL2.getConcatenatedTags(tempPos.maxVec2.midpoint, tempPos.maxVec2.range, tempPos.maxVec2.orientation);
			for(int i=0; i< profile.length; i++){
				profile[i] = profile[i]+tempaddToProfile.get(i);
			}
			int centerMidpoint = (temppair.BL1.equals(center) ? tempPos.maxVec1.midpoint : tempPos.maxVec2.midpoint);
			String centerorientation = (temppair.BL1.equals(center) ? tempPos.maxVec1.orientation : tempPos.maxVec2.orientation);
			int range = tempPos.maxVec1.range;
			for(BindingLocation jbl: seed){
				if(!jbl.equals(bl) && !jbl.equals(center)){
					CustomReturn cus = jbl.scanVecWithBl(center, centerMidpoint, centerorientation, range);
					tempaddToProfile = jbl.getConcatenatedTags(cus.maxvec.midpoint, range, cus.maxvec.orientation);
					for(int i=0; i< profile.length; i++){
						profile[i] = profile[i]+tempaddToProfile.get(i);
					}
				}
			}
			allseeds.put(bl, profile);
			
		}
		bestprofile = this.doPostAssessment(allseeds, conf);
		return bestprofile;
	}
	
	/**
	 * Much Simple implementation of Scheme 1. However, based on the assumption that the different offsets in the center don't mean much. 
	 * This scheme basically merges all the offsets that have the maximum similarity with any offset in the center. In other words, does not ensure that the offset in the 
	 * center is consistent with all comparisons. 
	 * Since only one profile is generated, does not call the post-assessment method.
	 * @return
	 */
	public int[] executeScheme2(Config conf){
		int[] profile=null;
		for(BindingLocation bl: this.seed){
			if(!bl.equals(center)){
				BLpair temppair = (allpairs.containsKey(new BLpair(center, bl)) ? new BLpair(center,bl): new BLpair(bl,center));
				CustomReturn tempPos = this.allpairs.get(temppair);
				List<Integer> tempaddToProfile = (temppair.BL1.equals(this.center) ? temppair.BL2.getConcatenatedTags(tempPos.maxVec2.midpoint, tempPos.maxVec2.range, tempPos.maxVec2.orientation): temppair.BL1.getConcatenatedTags(tempPos.maxVec1.midpoint, tempPos.maxVec1.range, tempPos.maxVec1.orientation));
				if(profile == null){
					profile = new int[tempaddToProfile.size()];
					for(int i=0; i<tempaddToProfile.size(); i++){
						profile[i] = tempaddToProfile.get(i);
					}
				}
				else{
					for(int i=0; i< profile.length; i++){
						profile[i] = profile[i]+tempaddToProfile.get(i);
					}
				}
		   }
		}
		//debug line
		//System.out.println(profile.length);
		CustomReturn centscan = center.scanConcVecWithBl(profile, conf.getIntSize());
		List<Integer> tempaddToProfile = this.center.getConcatenatedTags(centscan.maxvec.midpoint, centscan.maxvec.range, centscan.maxvec.orientation);
		for(int i=0; i<profile.length; i++){
			profile[i] = profile[i] + tempaddToProfile.get(i);
		}
		return profile;
	}
	
	public int[] executeScheme3(Config conf){
		int[] ret = new int[conf.getIntSize()*4];
		List<Node> workingTree = new LinkedList<Node>();
		
		List<BindingLocation> bllist = new ArrayList<BindingLocation>(this.pccpairwise.keySet());
		workingTree = this.getLinkedListFromListOfBl(bllist);
		
		double max_pcc = -10.0;
		
		do{
			max_pcc=-10.0;
			int left_pos =0;
			int right_pos=0;
			Node node_to_be_added=null;
			for(int i=0; i< workingTree.size()-1; i++){
				for(int j=i+1; j< workingTree.size(); j++){
					double temp_pcc;
					int temp_left_pos;
					int temp_right_pos;
					Node newnode;
					if(!workingTree.get(i).isleaf && !workingTree.get(j).isleaf){
						List<Integer> first = new ArrayList<Integer>();
						List<Integer> second =  new ArrayList<Integer>();
						for(int l=0; l< workingTree.get(i).composite.length; l++){
							first.add(workingTree.get(i).composite[l]);
							second.add(workingTree.get(j).composite[l]);
						}
						Pearson pcc_calculator = new Pearson(first, second);
						temp_pcc = pcc_calculator.doComparision();
						temp_left_pos = i;
						temp_right_pos = j;
						int[] newcomposite = new int[first.size()];
						for(int l=0; l< first.size(); l++){
							newcomposite[l] = first.get(l)+second.get(l);
						}
						newnode = new Node(workingTree.get(i),workingTree.get(j),newcomposite);
					}
					else if(workingTree.get(i).isleaf && workingTree.get(j).isleaf){
						CustomReturn cr = workingTree.get(i).leafbl.scanBlWithBl(workingTree.get(j).leafbl, conf.getIntSize());
						temp_pcc = cr.pcc;
						int[] newcomposite = new int[conf.getIntSize()*4];
						List<Integer> lefvec = workingTree.get(i).leafbl.getConcatenatedTags(cr.maxVec1.midpoint, cr.maxVec1.range, cr.maxVec1.orientation);
						List<Integer> rightvec = workingTree.get(j).leafbl.getConcatenatedTags(cr.maxVec2.midpoint, cr.maxVec2.range, cr.maxVec2.orientation);
						for(int l=0; l<lefvec.size(); l++){
							newcomposite[l] = lefvec.get(l)+rightvec.get(l);
						}
						temp_left_pos = i;
						temp_right_pos = j;
						newnode = new Node(workingTree.get(i), workingTree.get(j), newcomposite);
					}
					else{
						Node non_leaf_node = (workingTree.get(i).isleaf ? workingTree.get(j) : workingTree.get(i));
						Node leaf_node = (workingTree.get(i).isleaf ? workingTree.get(i) : workingTree.get(j));
						CustomReturn cr = leaf_node.leafbl.scanConcVecWithBl(non_leaf_node.composite, conf.getIntSize());
						temp_pcc = cr.pcc;
						int[] newcomposite = new int[conf.getIntSize()*4];
						List<Integer> vec = leaf_node.leafbl.getConcatenatedTags(cr.maxvec.midpoint, cr.maxvec.range, cr.maxvec.orientation);
						for(int l=0; l< vec.size(); l++){
							newcomposite[l] = non_leaf_node.composite[l]+vec.get(l);
						}
						temp_left_pos = i;
						temp_right_pos=j;
						newnode = new Node(workingTree.get(i), workingTree.get(j),newcomposite);
					}
					if(temp_pcc > max_pcc){
						node_to_be_added = newnode;
						left_pos = temp_left_pos;
						right_pos = temp_right_pos;
						max_pcc = temp_pcc;
					}
				}
			}
		
			if(max_pcc > conf.getSeedCutoff()){
				
				if(left_pos>right_pos){
					workingTree.remove(left_pos);
					workingTree.remove(right_pos);
				}
				else{
					workingTree.remove(right_pos);
					workingTree.remove(left_pos);
				}
				workingTree.add(node_to_be_added);
			}
			
			//debug line
			System.out.println(workingTree.size());
			System.out.println(max_pcc);
		} while(max_pcc>conf.getSeedCutoff() && workingTree.size()>0);
		
		int count_max=0;
		for(int l=0; l< workingTree.size(); l++){
			if(workingTree.get(l).count > count_max){
				count_max = workingTree.get(l).count;
				ret = workingTree.get(l).composite;
			}
			
		}
		System.out.println(count_max);
		
		//debug lines
		for(int m=0; m<ret.length; m++){
			System.out.println(m+"\t"+ret[m]);
		}
		
		
		return ret;
	}
	
	private List<Node> getLinkedListFromListOfBl(List<BindingLocation> bls){
		List<Node> ret = new LinkedList<Node>();
		for(int i=0; i< bls.size(); i++){
			Node tempnode = new Node(bls.get(i));
			ret.add(tempnode);
		}
		return ret;
	}
	

	
	//Accessors
	public int getNoInSeed(){ return seed.size();}
	public int getNoRemoved(){return removed.size();}

}
