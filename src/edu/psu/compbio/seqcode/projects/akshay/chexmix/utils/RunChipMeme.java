package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;

public class RunChipMeme {
	public List<BindingLocation> locations;
	
	public RunChipMeme(List<BindingLocation> locations) {
		this.locations = locations;
	}
	
	public void printFasta(){
		for(int i=0; i<locations.size(); i++){
			System.out.println("loc"+i+"\n");
			System.out.println(locations.get(i).seqpos+"\n");
			
		}
	}
	
	public static void main(String[] args) throws IOException{
		BufferedReader br = null;
		br =  new BufferedReader(new FileReader("/gpfs/home/auk262/scratch/test_input"));
		RunChipMeme chipdriver = null;
		List<BindingLocation> locationlist = new ArrayList<BindingLocation>();
		
		String currentline = br.readLine();
		while(currentline != null){
			String[] pieces = currentline.split("\t");
			int tempMidpoint = (Integer.parseInt(pieces[4])-Integer.parseInt(pieces[3]))/2;
			String tempChr = pieces[0];
			BindingLocation temploc = new BindingLocation(tempMidpoint, tempChr, 60);
			locationlist.add(temploc);
			currentline = br.readLine();
		}
		br.close();
		chipdriver = new RunChipMeme(locationlist);
		chipdriver.printFasta();
		
	}

}
