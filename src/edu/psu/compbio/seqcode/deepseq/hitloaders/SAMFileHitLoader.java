package edu.psu.compbio.seqcode.deepseq.hitloaders;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import edu.psu.compbio.seqcode.deepseq.HitPair;
import edu.psu.compbio.seqcode.deepseq.Read;
import edu.psu.compbio.seqcode.deepseq.ReadHit;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;


/**
 * SAMFileHitLoader: A FileHitLoader for SAM and BAM files.
 * Accounts for uniqueness of hits according to user-specified option. 
 * Ignores secondary & supplementary (i.e. chimeric) alignments.
 * @author mahony
 *
 */
public class SAMFileHitLoader extends FileHitLoader{

	private boolean useChimericReads=false; //Ignore chimeric mappings for now. 
	
	public SAMFileHitLoader(File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadRead2, boolean loadPairs) {
    	super(f, nonUnique, true, false, loadRead2, loadPairs);
    	if(!loadT1Reads || loadT2Reads)
			System.err.println("SAMFileHitLoader: You asked to load only Type1 or Type2 reads, we do not yet load this information from SAM format.");
    }
    
    /**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 * Loads pairs to hitPairsList
	 */
	public void sourceAllHits() {
		this.initialize();
		SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		SamReader reader = factory.open(file);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
		String lastread = null;
		while (iter.hasNext()) {
		    SAMRecord record = iter.next();
		    
		    if(record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary() && !useChimericReads){continue;}
		    if(record.getReadPairedFlag() && !record.getFirstOfPairFlag() && !loadRead2){continue;}
		    	
		    if (lastread == null || !lastread.equals(record.getReadName())) {
		    	processRead(byRead);
		    	byRead.clear();
		    }
		    lastread = record.getReadName();
		    
		    byRead.add(record); //Filter by first or second of pair here if loading by type1/2?

		    //load pair if this is a first mate, congruent, proper pair
		    if(loadPairs && record.getFirstOfPairFlag() && record.getProperPairFlag()){
		    	boolean neg = record.getReadNegativeStrandFlag();
                boolean mateneg = record.getMateNegativeStrandFlag();
                HitPair hp = new HitPair((neg ? record.getAlignmentEnd() : record.getAlignmentStart()),
                		record.getMateReferenceName().replaceFirst("^chr", ""),
                		(mateneg ? record.getMateAlignmentStart()+record.getReadLength()-1 : record.getMateAlignmentStart()), 
                		mateneg ? 1 : 0,
                		1);
                addPair(record.getReferenceName().replaceFirst("^chr", ""), neg ? '-':'+', hp);
		    }
		}
		processRead(byRead);
		iter.close();
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }//end of sourceAllHits method
    
    protected void processRead(Collection<SAMRecord> records) {
        int mapcount = records.size();
        if(mapcount == 0)
            return;
        if(!useNonUnique && mapcount > 1)
            return;
        
        float weight = 1 / ((float)mapcount);
        Read currRead = new Read();
		for (SAMRecord record : records) {
		    int start = record.getAlignmentStart();
		    int end = record.getAlignmentEnd();
		    ReadHit currHit = new ReadHit(
						  record.getReferenceName().replaceFirst("^chr", ""), 
						  start, end, 
						  record.getReadNegativeStrandFlag() ? '-' : '+',
						  weight);
		    currRead.addHit(currHit);
		}currRead.setNumHits(mapcount);
		addHits(currRead);
    }//end of processRead
}
