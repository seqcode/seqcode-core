package edu.psu.compbio.seqcode.deepseq.hitloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import edu.psu.compbio.seqcode.deepseq.Read;
import edu.psu.compbio.seqcode.deepseq.ReadHit;

/**
 * IDXFileHitLoader: a FileHitLoader from Pugh lab IDX format
 * 
 * @author mahony
 *
 */
public class IDXFileHitLoader extends FileHitLoader {

	public IDXFileHitLoader(File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadPairs){
		super(f, nonUnique, true, false, false, loadPairs);
		if(!loadT1Reads || loadT2Reads)
			System.err.println("IDXFileHitLoader: You asked to load only Type1 or Type2 reads, but IDX cannot discriminate between reads for single-end hits.");
		if(loadPairs)
			System.err.println("IDXFileHitLoader: You asked to load pairs, but IDX cannot represent paired read data.");
	}
	
	/**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 * Nothing stored to hitPairsList, since IDX cannot store pairs. 
	 */
	public void sourceAllHits() {
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
            			chr=tmp[0].replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", "");
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
	
}//end of IDXFileHitLoader class
