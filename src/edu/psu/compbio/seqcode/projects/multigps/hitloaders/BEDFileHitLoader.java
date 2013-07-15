package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import edu.psu.compbio.seqcode.projects.multigps.framework.*;

public class BEDFileHitLoader extends FileHitLoader {

	public BEDFileHitLoader(File f, boolean nonUnique) {
    	super(f, nonUnique);
	}
	
	/**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 */
	public void sourceReads() {
		this.initialize();
		try {
			totalHits=0;
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line;
			float currReadHitCount=0;
			Read currRead=null;
	        while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\s+");
		            String chr="."; char strand = '.';
		            int start=0, end=0;
		           
		            if(currRead!=null){
		            	currRead.setNumHits(currReadHitCount);
	            		//Add the hits to the data structure
	            		addHits(currRead);
	           			currRead=null;
	           		}
	            	currReadHitCount=1;            			
	            	try{
            			chr = words[0];
            			String[] tmp = chr.split("\\.");
            			chr=tmp[0].replaceFirst("chr", "");
            			chr=chr.replaceFirst("^>", "");
	// http://genome.ucsc.edu/FAQ/FAQformat.html#format1
	//BED format is half open - The chromEnd base is not included  
	// For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    //BED format is also 0-based, and we want 1-based
	            		start = new Integer(words[1]).intValue()+1;
	            		end = new Integer(words[2]).intValue();
	           			strand = words[5].charAt(0);
	    				ReadHit currHit = new ReadHit(chr, start, end, strand);
	    				currRead = new Read();
	    				currRead.addHit(currHit);
	            	} catch (NumberFormatException e){
	            		// skip reading this line for header or comment lines
	           		}
	    		}
            }
	        if(currRead!=null){
    			currRead.setNumHits(currReadHitCount);
    			//Add the hits to the data structure
    			addHits(currRead);
    		}
	        reader.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of BEDFileReader class
