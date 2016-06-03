package org.seqcode.deepseq.hitloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.seqcode.deepseq.Read;
import org.seqcode.deepseq.ReadHit;


/**
 * BowtieFileHitLoader: a FileHitLoader from Bowtie native format
 * 
 * @author mahony
 *
 */
public class BowtieFileHitLoader extends FileHitLoader {

	public BowtieFileHitLoader(File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadPairs){
		super(f, nonUnique, true, false, false, loadPairs);
		if(!loadT1Reads || loadT2Reads)
			System.err.println("BowtieFileHitLoader: You asked to load only Type1 or Type2 reads, we do not load this information from Bowtie native format.");
		if(loadPairs)
			System.err.println("BowtieFileHitLoader: You asked to load pairs, but we do not load paired read data from Bowtie native format.");
	}	
		
	/**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 * Nothing loaded to hitPairsList, since we do not load this type of data from Bowtie native
	 */
	public void sourceAllHits() {
		this.initialize();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line, lastID="";
			double currReadHitCount=0;
			Read currRead=null;
			int readLength=-1;
	        while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if (line.length()==0) continue;
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            String chr="."; char strand = '.';
		            int start=0, end=0;
		            
		            String ID = words[0];
		            if(readLength==-1)
	    				readLength = words[4].length();
		            
		            boolean newRead=true;
	            	if(ID.equals(lastID)){
	            		currReadHitCount++;
	            		newRead=false;
	            	}else{
	            		if(currRead!=null){
	            			if(currRead.getNumHits()==1 || useNonUnique){
	        	        		//Add the hits to the data structure
	        	        		addHits(currRead);
	        	        	}currRead=null;
	            		}
	            		currReadHitCount=1;            			
	            	}
	            	int mis=0;
	            	if(words.length>7 && words[7].length()>1){
	            		mis = words[7].split(",").length;
	            	}
	            	chr = words[2];
        			String[] tmp = chr.split("\\.");
        			chr=tmp[0].replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", "");
        			chr=chr.replaceFirst("^>", "");
        			start = new Integer(words[3]).intValue()+1; //Bowtie raw output is 0-based, we want 1-based
        			end =start+readLength-1;
        			strand = words[1].charAt(0);
					ReadHit currHit = new ReadHit(chr, start, end, strand, 1);
					if(newRead || currRead==null){
						currRead = new Read();
	            	}
					currRead.addHit(currHit);
    			
	            	lastID=ID;
	        	}
            }
	        if(currRead!=null){
	        	if(currRead.getNumHits()==1 || useNonUnique){
	        		//Add the hits to the data structure
	        		addHits(currRead);
	        	}
    		}
	        reader.close();
	 
		} catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of BowtieFileReader class
