package edu.psu.compbio.seqcode.gse.projects.readdb;

import java.io.*;
import java.util.*;

import org.apache.commons.cli.*;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

/**
 * Reads SAM or BAM data on stdin.
 * Produces a file on stdout in the format expected by ImportHits.
 * 
 * Options:	--uniquehits (flag to only print 1:1 read to hit mappings)
 * 			--pairedend (flag to print pairs)
 * 			--junctions (flag to print junction mapping reads as pairs)
 * 
 */

public class Bowtie2SAMToReadDB {

    public static boolean uniqueOnly;
    public static boolean inclPairedEnd;
    public static boolean inclJunction;
    
    public static boolean lastFirstMateUnique=false; 

    public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("u","uniquehits",false,"only output hits with a single mapping");
        options.addOption("p","pairedend",false,"output paired-end hits");
        options.addOption("j","junctions",false,"output junction mapping reads (reads with a single gap)");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );            
    	uniqueOnly = cl.hasOption("uniquehits");
    	inclPairedEnd = cl.hasOption("pairedend");
    	inclJunction = cl.hasOption("junctions");
    	SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        SAMFileReader reader = new SAMFileReader(System.in);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            if (record.getReadUnmappedFlag()) {continue; }
            processRecord(record);
        }
        iter.close();
        reader.close();
    }       
    public static void processRecord(SAMRecord record) {

    	int primAScore = record.getIntegerAttribute("AS");
    	int secAScore=-1000000;
    	if(record.getIntegerAttribute("XS")!=null)
    		secAScore = record.getIntegerAttribute("XS");
    	
    	if(inclPairedEnd || inclJunction){
    		boolean currUnique = primAScore > secAScore ? true : false;
        	float weight = 1;  //Fix this if using bowtie2 to produce multiple mappings for each read
    	    
    		/*
    		 * Only accept proper, congruent pairs.
    		 * It also assumes that the left and right mates have the same length, 
    		 * and that there are no gaps in the second mate alignment (SAM doesn't store the paired read's end)
    		 * Note: if you change this, you may have to change the SAMStats output also
    		 */
        	if(inclPairedEnd){
        		if(record.getFirstOfPairFlag() && record.getProperPairFlag() && record.getStringAttribute("YT").equals("CP")){
        			lastFirstMateUnique = currUnique;
        		}else if(record.getSecondOfPairFlag() && record.getProperPairFlag() && record.getStringAttribute("YT").equals("CP")){
    			
        			if(!uniqueOnly || (currUnique || lastFirstMateUnique)){
		    			//Print
		                boolean neg = record.getReadNegativeStrandFlag();
		                boolean mateneg = record.getMateNegativeStrandFlag();
		                String len = record.getReadLength() + "\t";
		                System.out.println(
		                		record.getMateReferenceName() + "\t" +
		                		(mateneg ? 
		                			record.getMateAlignmentStart()+record.getReadLength()-1 : 
		                			record.getMateAlignmentStart()) + "\t" +
		                		(mateneg ? "-\t" : "+\t") +
		                		len +
		                                
		                        record.getReferenceName() + "\t" +
		                       	(neg ? 
		                       		record.getAlignmentEnd() : 
		                       		record.getAlignmentStart()) + "\t" +
		                        (neg ? "-\t" : "+\t") + 
		                        len +
		                        
		                        weight +"\t"+
		                        1);
        			}
    			}
        	}
			/*
    		 * Outputs as paired alignments those reads that are aligned in >2 blocks
    		 * Note: if you change this, you may have to change the SAMStats output also
    		 */
    		if(inclJunction){
    			if(!uniqueOnly || currUnique){
		    		List<AlignmentBlock> blocks = record.getAlignmentBlocks();
		    		if(blocks.size()>=2){
		    			for(int ab=0; ab<blocks.size()-1; ab++){
			    			AlignmentBlock lBlock = blocks.get(ab);
			    		   	int lStart = lBlock.getReferenceStart();
			    		   	int lEnd = lStart + lBlock.getLength()-1;
			    		   	int lLen = lBlock.getLength();
			    		   	AlignmentBlock rBlock = blocks.get(ab+1);
			    		   	int rStart = rBlock.getReferenceStart();
			    		   	int rEnd = rStart + rBlock.getLength()-1;
			    		   	int rLen = rBlock.getLength();
			                boolean neg = record.getReadNegativeStrandFlag();
			                String refname = record.getReferenceName() + "\t";
			    		   	System.out.println(
			                                   refname +
			                                   (neg ? lEnd : lStart) + "\t" +
			                                   (neg ? "-\t" : "+\t") +
			                                   lLen + "\t" +
			                                   refname + 
			                                   (neg ? rEnd : rStart) + "\t" +
			                                   (neg ? "-\t" : "+\t") +
			                                   rLen + "\t" +
			                                   weight +"\t"+
			                                   0);
		    			}
		    		}
    			}
    		}
    	}else{ //Just output reads (ignore alignment blocks for now)
    		
        	if (uniqueOnly && primAScore == secAScore) {
                return;
            }
        	float weight = 1;  //Fix this if using bowtie2 to produce multiple mappings for each read
    	    
    		System.out.println(String.format("%s\t%d\t%s\t%d\t%f",
                    record.getReferenceName(),
                    record.getReadNegativeStrandFlag() ? 
                    record.getAlignmentEnd() : 
                    record.getAlignmentStart(),
                    record.getReadNegativeStrandFlag() ? "-" : "+",
                    record.getReadLength(),
                    weight));
    	}                                 
    }
}