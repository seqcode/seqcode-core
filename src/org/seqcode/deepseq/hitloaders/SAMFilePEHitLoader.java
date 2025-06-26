package org.seqcode.deepseq.hitloaders;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Collection;

import org.seqcode.deepseq.HitPair;
import org.seqcode.deepseq.Read;
import org.seqcode.deepseq.ReadHit;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;


/**
 * SAMFilePairHitLoader: A FileHitLoader for pair-end SAM and BAM files.
 * Accounts for uniqueness of hits according to user-specified option. 
 * Ignores secondary & supplementary (i.e. chimeric) alignments.
 * Assume SAM/BAM file is sorted
 * 
 * @author jianyu
 *
 */
public class SAMFilePEHitLoader extends FileHitLoader{

	private boolean useChimericReads=false; //Ignore chimeric mappings for now. 
	
	public SAMFilePEHitLoader(File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads) {
    	super(f, nonUnique, true, false, true, true);
    	if(!loadT1Reads || loadT2Reads)
			System.err.println("SAMFileHitLoader: You asked to load only Type1 or Type2 reads, we do not yet load this information from SAM format.");
    }
    
    /**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads pairs to hitPairsList
	 */
	public void sourceAllHits() {
		this.initialize();
		SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
		            		  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		SamReader reader = factory.open(file);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		Map<String, SAMRecord> byRead = new HashMap<String, SAMRecord>();
		
		SAMRecord r1Record;
		SAMRecord r2Record;
		while (iter.hasNext()) {
		    SAMRecord record = iter.next();
		    
		    if(record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary() && !useChimericReads){continue;}
		    
		    if(byRead.containsKey(record.getReadName())) {
		    	if (record.getFirstOfPairFlag()) {
		    		r1Record = record;
		    		r2Record = byRead.remove(record.getReadName());
		    	} else {
		    		r1Record = byRead.remove(record.getReadName());
		    		r2Record = record;
		    	}
		    	boolean neg = r1Record.getReadNegativeStrandFlag();
		    	boolean mateneg = r2Record.getReadNegativeStrandFlag();
		    	
		    	// skip the pair if it's not aligned correctly
		    	if (!r1Record.getReferenceName().equals(r2Record.getReferenceName())) {continue;}
		    	if (neg == mateneg) {continue;}

		    	HitPair hp = new HitPair((neg ? r1Record.getAlignmentEnd() : r1Record.getAlignmentStart()),
		    			r2Record.getReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""),
                		(mateneg ? r2Record.getAlignmentEnd() : r2Record.getAlignmentStart()), 
                		mateneg ? 1 : 0,
                		1);
		    	addPair(r1Record.getReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""), neg ? '-':'+', hp);
		    			    	
		    } else {
		    	byRead.put(record.getReadName(), record);
		    }
		    // add the hit
//		    ReadHit currHit = new ReadHit(
//						  record.getReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""), 
//						  record.getAlignmentStart(), record.getAlignmentEnd(), 
//						  record.getReadNegativeStrandFlag() ? '-' : '+',
//						  1);
//		    Read currRead = new Read();
//		    currRead.addHit(currHit);
//		    currRead.setNumHits(1);
//		    addHits(currRead);
	    }

		iter.close();
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }//end of sourceAllHits method
  
}
