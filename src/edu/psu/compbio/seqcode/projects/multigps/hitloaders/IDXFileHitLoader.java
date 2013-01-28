package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import edu.psu.compbio.seqcode.projects.multigps.framework.*;

public class IDXFileHitLoader extends FileHitLoader {

	public IDXFileHitLoader(File f, boolean nonUnique) {
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
			Read posRead=null, negRead=null;
	        while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	String[] words = line.split("\\s+");
	        	if(line.charAt(0)!='#' && words.length>=4 && !words[0].equals("chrom")){
		            String chr="."; char strand = '.';
		            int fivePrime=0;
		            float posWeight=0, negWeight=0;
		            try{
            			chr = words[0];
            			String[] tmp = chr.split("\\.");
            			chr=tmp[0].replaceFirst("chr", "");
            			chr=chr.replaceFirst("^>", "");
            			fivePrime = new Integer(words[1]).intValue();
            			posWeight = new Integer(words[2]).floatValue();
            			negWeight = new Integer(words[3]).floatValue();
            			if(posWeight>0){
            				ReadHit posHit = new ReadHit(chr, fivePrime, fivePrime, '+', posWeight);
            				posRead = new Read();
            				//Add the hit, but don't update weight since the weight is used here for representing the read count
            				posRead.addHit(posHit, false);
            				addHits(posRead);
            			}
            			if(negWeight>0){
            				ReadHit negHit = new ReadHit(chr, fivePrime, fivePrime, '-', negWeight);
            				negRead = new Read();
            				//Add the hit, but don't update weight since the weight is used here for representing the read count
            				negRead.addHit(negHit, false);
            				addHits(negRead);
            			}
	            	} catch (NumberFormatException e){
	            		// skip reading this line for header or comment lines
	           		}
	    		}
            }
	        reader.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of BEDFileReader class
