package edu.psu.compbio.seqcode.projects.akshay.utils;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class ExtractPairedMidPoints {
	
	
	
	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		for(String s : ap.getKeys()){
			if(s.startsWith("rdb")){
				String condrep = s.replaceFirst("rdbexpt", "");
    			
			}
		}
		
		
		
	}
	

}
