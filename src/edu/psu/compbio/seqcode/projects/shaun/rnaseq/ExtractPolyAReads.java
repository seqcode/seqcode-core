package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import java.io.File;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

public class ExtractPolyAReads{
	File inFile=null;
	int minARun=5;
	int minReadLength = 20;
	String minARunString="";
	String minTRunString="";
	
    public ExtractPolyAReads(File f) {
    	inFile = f;
    	for(int a=0; a<minARun; a++){
    		minARunString = minARunString.concat("A");
    		minTRunString = minTRunString.concat("T");
    	}
    }

    /**
	 * Find:
	 * 	 - Unmapped reads
	 *   - That end in a run of As or begin with a run of Ts
	 */
	public void execute() {
		SAMFileReader reader = new SAMFileReader(inFile);
		
		reader.setValidationStringency(ValidationStringency.SILENT);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
		    SAMRecord record = iter.next();
		    
		    if (record.getReadUnmappedFlag()) {
		    	
		    	String seq = record.getReadString();
		    	if(seq.endsWith(minARunString) || seq.startsWith(minTRunString)){
		    		String qual = record.getBaseQualityString();
		    		String newSeq="";
		    		String newQual="";
		    		if(seq.endsWith(minARunString)){
		    			//Iterate to the start of the non-polyA
		    			int x;
		    			for(x=seq.length(); x>0 && seq.charAt(x)=='A'; x--){}
		    			
		    			newSeq = seq.substring(0, x+1);
		    			newQual = qual.substring(0, x+1);
		    		}else{
		    			//Iterate to the start of the non-polyA
		    			int x;
		    			for(x=0; x<seq.length() && seq.charAt(x)=='T'; x++){}

					//For poly-T starting reads, return the reverse complement seq and the reverse qual
                                        newSeq = SequenceUtils.reverseComplement(seq.substring(x, seq.length()));

                                        String qualSeg = qual.substring(x, seq.length());
                                        StringBuilder sb = new StringBuilder();
                                        for(int q = qualSeg.length()-1; q>= 0; q--)
                                            sb.append(qualSeg.charAt(q));
                                        newQual = sb.toString();
		    		}
		    		//System.err.println("Old:\t"+seq+"\nNew:\t"+newSeq);
		    		
		    		if(newSeq.length()>=minReadLength)
		    			System.out.println("@"+record.getReadName()+"\n"+newSeq+"\n+"+record.getReadName()+"\n"+newQual);
		    	}
			}	
		}
		iter.close();
		reader.close();
    }
	
	public static void main(String[] args){
		if(args.length==0){
			System.err.println("ExtractPolyAReads:\n" +
					"\t--in [SAM/BAM file]\n");
		}else{
			ArgParser ap = new ArgParser(args);
			
			if(ap.hasKey("in")){
				String inFile = ap.getKeyValue("in");
				
				ExtractPolyAReads extractor = new ExtractPolyAReads(new File(inFile));
				extractor.execute();
			}else{
				System.err.println("No input file specified");
			}
		}
	}
}
