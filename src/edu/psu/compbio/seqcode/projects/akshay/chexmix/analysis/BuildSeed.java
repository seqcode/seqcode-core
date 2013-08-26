package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;


public class BuildSeed {
	
	public Map<BLpair,CustomReturn> allpairs = new HashMap<BLpair, CustomReturn>();
	public Map<BindingLocation,Integer> removed = new HashMap<BindingLocation,Integer>();
	public Map<BindingLocation,Double> pccpairwise = new HashMap<BindingLocation,Double>();
	public double cutoff;
	public int smoothsize;
	public BuildSeed(List<BindingLocation> topBlList, int range, int smoothsize, double cutoff) {
		this.cutoff = cutoff;
		this.smoothsize = smoothsize;
		for(int i=0; i<topBlList.size(); i++){
			for(int j=0; j< topBlList.size(); j++ ){
				if(i != j){
					BLpair pair = new BLpair(topBlList.get(i),topBlList.get(j));
					CustomReturn pccmax = topBlList.get(i).scanTwoBLs(topBlList.get(j), range, smoothsize);
					allpairs.put(pair, pccmax);
				}
			}
		}
	}
	
	public BindingLocation getCenter(){
		for (BLpair pair: allpairs.keySet()){
			if(!removed.containsKey(pair.BL1) && !removed.containsKey(pair.BL2)){
				if(pccpairwise.containsKey(pair.BL1)){
					pccpairwise.put(pair.BL1, pccpairwise.get(pair.BL1) + allpairs.get(pair).pcc);
				}
				else{
					pccpairwise.put(pair.BL1, allpairs.get(pair).pcc);
				}
			}
		}
		
		BindingLocation ret=null;
		double max=-1000;
		
		for(BindingLocation bl : pccpairwise.keySet()){
			if(!removed.containsKey(bl)){
				if(pccpairwise.get(bl)>max){
					max =pccpairwise.get(bl);
					ret = bl;
				}
			}
		}
		return ret;
	}
	
	public void updateRemoved(BindingLocation center){
		for(BLpair pair : allpairs.keySet()){
			if(pair.BL1 == center){
				if(allpairs.get(pair).pcc < this.cutoff){
					removed.put(pair.BL2, 1);
				}
			}
		}
		
	}
	
	public void printSeed(){
		BindingLocation center = this.getCenter();
		List<Integer> cumtags=null;
		int count=0;
		for(BLpair pair: allpairs.keySet()){
			
			
			if(pair.BL1 == center && !removed.containsKey(pair.BL2)){
				Vec max= allpairs.get(pair).maxVec2;
				List<Integer> tags = pair.BL2.getConcatenatedTags(max.midpoint, max.chr, max.range, max.orientation, this.smoothsize);
				if(cumtags == null){
					cumtags = tags;
				}
				else{
					for(int i=0; i<tags.size(); i++){
					cumtags.add(cumtags.get(i)+tags.get(i));
					}
				}
			}
		}
		System.out.println("No of location in seed:"+count);
		for(int i:cumtags){
			System.out.println(i);
		}
		
	}
	/**
	 * this main body is only to test the buildseed class.
	 * @param args
	 * @throws IOException
	 */
	
	public static void main(String[] args) throws IOException{
		System.out.println("Loading the hitloader");
		LoadTags loader = new LoadTags();
		loader.loadHits("/gpfs/home/auk262/scratchFoxA2_07-633_liver_-_-_-_XO_kaz1-S001_Pugh4020mm10.idx", "idx", false);
		BufferedReader br = null;
		File tempdir = new File("temp");
		tempdir.mkdir();
		br =  new BufferedReader(new FileReader("/gpfs/home/auk262/scratch/test_input"));
		List<BindingLocation> locationlist = new ArrayList<BindingLocation>();
		String currentline = br.readLine();
		while(currentline != null){
			String[] pieces = currentline.split("\t");
			int tempMidpoint = Integer.parseInt(pieces[3])+((Integer.parseInt(pieces[4])-Integer.parseInt(pieces[3]))/2);
			String tempChr = pieces[0];
			System.out.println("Creating Binding Location");
			BindingLocation temploc = new BindingLocation(tempMidpoint, tempChr, 300);
			System.out.println("Filling Tags");
			temploc.filltags(loader);
			locationlist.add(temploc);
			currentline = br.readLine();
		}
		System.out.println(locationlist.get(1).vecpos.tags);
		br.close();
		
		BuildSeed driver = new BuildSeed(locationlist, 60, 0, 0.5);
		for (int i=0; i<5; i++){
			System.out.println("Getting Center");
			BindingLocation tempcenter = driver.getCenter();
			driver.updateRemoved(tempcenter);
		}
		driver.printSeed();
		
	}
	
}
