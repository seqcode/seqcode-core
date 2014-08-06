package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class ExtractPolyAReadsFromFastq{
	File fqFile=null;
	File samFile=null;
	boolean mappedFilter=false;
	int minARun=5;
	int minReadLength = 20;
	String minARunString="";
	String minTRunString="";
	boolean leftReads=true;
	
    public ExtractPolyAReadsFromFastq(File fastq, File sam, boolean left) {
    	fqFile = fastq;
    	samFile = sam;
    	leftReads = left;
    	for(int a=0; a<minARun; a++){
    		minARunString = minARunString.concat("A");
    		minTRunString = minTRunString.concat("T");
    	}
    	if(samFile!=null)
    		mappedFilter=true;
    }

    /**
	 * Find:
	 * 	 - Unmapped reads
	 *   - That end in a run of As or begin with a run of Ts
	 */
	public void execute() {
		HashMap<String, Integer> mappedReads=null;
		if(mappedFilter){
			mappedReads =  new HashMap<String,Integer>(150000000);
			SamReaderFactory factory =
			          SamReaderFactory.makeDefault()
			              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			              .validationStringency(ValidationStringency.SILENT);
			SamReader sreader = factory.open(samFile);
			CloseableIterator<SAMRecord> iter = sreader.iterator();
			while (iter.hasNext()) {
			    SAMRecord record = iter.next();
			    
			    if (!record.getReadUnmappedFlag()) {
			    	if((leftReads && record.getFirstOfPairFlag()) || (!leftReads && record.getSecondOfPairFlag())){
			    		String name = record.getReadName();
			    		if(leftReads)
			    			name = name.concat("/1");
			    		else
			    			name = name.concat("/2");
			    		mappedReads.put(name, 1);
			    	}
			    }
			}
			iter.close();
			try {
				sreader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
			System.err.println(mappedReads.keySet().size()+" mapped reads in the SAM file.");
		}
		int polyACount =0;
		int totalReads =0;
		int unmapped=0;
		//Iterate through FASTQ file
		try{
		    if(!fqFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(fqFile));
			String line;
			while ((line = reader.readLine()) != null) {
	        	totalReads++;
			line = line.trim();
			
			String name = line.replaceFirst("@", "");
			String seq = reader.readLine().trim();
			String name2 = reader.readLine().trim();
			String qual = reader.readLine().trim();
			
			if(!mappedFilter || !mappedReads.containsKey(name)){
			    unmapped++;
			    if(seq.endsWith(minARunString) || seq.startsWith(minTRunString)){
	    	    		String newSeq="";
	    	    		String newQual="";
	    	    		if(seq.endsWith(minARunString)){
	    	    			//Iterate to the start of the non-polyA
	    	    			int x;
	    	    			for(x=seq.length()-1; x>0 && seq.charAt(x)=='A'; x--){}
	    	    			
	    	    			newSeq = seq.substring(0, x+1);
	    	    			newQual = qual.substring(0, x+1);
	    	    		}else{
	    	    			//Iterate to the start of the non-polyA
	    	    			int x;
	    	    			for(x=0; x<seq.length()-1 && seq.charAt(x)=='T'; x++){}
	    	    			
	    	    			//For poly-T starting reads, return the reverse complement seq and the reverse qual
	    	    			newSeq = SequenceUtils.reverseComplement(seq.substring(x, seq.length()));
					
	    	    			String qualSeg = qual.substring(x, seq.length());
	    	    			StringBuilder sb = new StringBuilder();
	    	    			for(int q = qualSeg.length()-1; q>= 0; q--)
	    	    				sb.append(qualSeg.charAt(q));
	    	    				newQual = sb.toString();
	    	    		}
	    	    			//System.err.println("Old:\t"+seq+"\t"+qual+"\nNew:\t"+newSeq+"\t"+newQual);
	    	    		
	    	    		if(newSeq.length()>=minReadLength){
	    	    			System.out.println("@"+name+"\n"+newSeq+"\n"+name2+"\n"+newQual);
	    	    			polyACount++;
	    	    		}
	    	    	}
	            }
	        }
	        System.err.println(totalReads+" total reads in the FASTQ file.");
	        System.err.println(unmapped+" unmapped reads.");
	        System.err.println(polyACount+" potential polyA reads.");
	        reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}	
	
	public static void main(String[] args){
		if(args.length==0){
			System.err.println("ExtractPolyAReadsFromFastq:\tUse this version with SAM files that don't include the unmapped reads.\n" +
					"\t--sam [SAM/BAM file]\n" +
					"\t--fq [fastq file]\n" +
					"\t--side <left/right>\n");
		}else{
			ArgParser ap = new ArgParser(args);
			
			boolean leftReads =true;
			if(ap.hasKey("fq")){
				String samFile = null;
				if(ap.hasKey("sam"))
					ap.getKeyValue("sam");
				String fqFile = ap.getKeyValue("fq");
				if(ap.hasKey("side"))
					leftReads = ap.getKeyValue("side").equals("left");
				
				ExtractPolyAReadsFromFastq extractor = new ExtractPolyAReadsFromFastq(new File(fqFile), samFile==null ? null : new File(samFile), leftReads);
				extractor.execute();
			}else{
				System.err.println("No input file specified");
			}
		}
	}
}
