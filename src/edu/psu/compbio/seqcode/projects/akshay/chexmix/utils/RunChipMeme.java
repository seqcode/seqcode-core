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
			System.out.println("loc"+i);
			System.out.println(locations.get(i).seqpos.sequence);
			
		}
	}
	
	

}
