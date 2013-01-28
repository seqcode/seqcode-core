package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import edu.psu.compbio.seqcode.projects.multigps.framework.*;

public class NovoFileHitLoader extends FileHitLoader {

	public NovoFileHitLoader(File f, boolean nonUnique) {
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
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            String chr="."; char strand = '.';
		            int start=0, end=0;
		            
		            String ID = words[0];
		            if(readLength==-1)
	    				readLength = words[2].length();
		            
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
	            	String tag = words[4];
	            	if(tag.equals("U") || (useNonUnique && words.length>9 && tag.charAt(0)=='R')){
	            		int mis=0;
		            	if(words.length>13 && words[13].length()>1){
		            		mis = words[7].split(" ").length;
		            	}
            			chr = words[7];
            			String[] tmp = chr.split("\\.");
            			chr=tmp[0].replaceFirst("chr", "");
            			chr=chr.replaceFirst("^>", "");
            			start = new Integer(words[8]).intValue();
            			end =start+readLength-1;
            			strand = words[9].equals("F") ? '+' : '-';
    					ReadHit currHit = new ReadHit(chr, start, end, strand, 1);
    					if(!ID.equals(lastID) || currRead==null){
    						currRead = new Read();
    	            	}currRead.addHit(currHit);
	    			}
	            	lastID=ID;
	        	}
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

}//end of NovoFileReader class
