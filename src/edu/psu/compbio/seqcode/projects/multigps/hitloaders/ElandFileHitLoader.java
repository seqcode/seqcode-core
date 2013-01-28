package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import edu.psu.compbio.seqcode.projects.multigps.framework.*;

public class ElandFileHitLoader extends FileHitLoader {

	public ElandFileHitLoader(File f, boolean nonUnique) {
    	super(f, nonUnique);
	}	
		
	/**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 */
	public void sourceReads() {
		this.initialize();
		try {
			int readLength=-1;
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line, lastID="";
			float currReadHitCount=0;
			Read currRead=null;
	        while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	            String[] words = line.split("\\t");
	            String chr="."; char strand = '.';
	            int start=0, end=0;
	            
	            String ID = words[0];
	            if(readLength==-1)
    				readLength = words[1].length();
	            
            	if(ID.equals(lastID)){
            		currReadHitCount++;
            	}else{
            		if(currRead!=null){
            			currRead.setNumHits(currReadHitCount);
            			//Add the hits to the data structure
            			addHits(currRead);
            			currRead=null;
            		}
            		currReadHitCount=1;            			
            	}
            	String tag = words[2];
            	if(tag.charAt(0)=='U' || (useNonUnique && words.length>8 && tag.charAt(0)=='R')){
            		int mis = Integer.parseInt(tag.substring(1, 2));
        			chr = words[6];
        			String[] tmp = chr.split("\\.");
        			chr=tmp[0].replaceFirst("^chr", "");
        			start = new Integer(words[7]).intValue();
        			end =start+readLength-1;
        			strand = words[8].equals("F") ? '+' : '-';
					ReadHit currHit = new ReadHit(chr, start, end, strand, 1);
					if(!ID.equals(lastID) || currRead==null){
						currRead = new Read();
	            	}currRead.addHit(currHit);
    			}
            	lastID=ID;
            }
	        if(currRead!=null){
    			currRead.setNumHits(currReadHitCount);
    			//Add the hits to the data structure
    			addHits(currRead);
    		}
	        reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of ElandFileReader class
