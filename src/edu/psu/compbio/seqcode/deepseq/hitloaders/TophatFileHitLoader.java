package edu.psu.compbio.seqcode.deepseq.hitloaders;


import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.HitPair;
import edu.psu.compbio.seqcode.deepseq.Read;
import edu.psu.compbio.seqcode.deepseq.ReadHit;

/**
 * TophatFileHitLoader: A FileHitLoader for Tophat SAM/BAM output.
 * Each alignment block is loaded as a separate hit. Ignores uniqueness. 
 * @author mahony
 *
 */
public class TophatFileHitLoader extends FileHitLoader{

    public TophatFileHitLoader(File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadRead2, boolean loadPairs) {
    	super(f, nonUnique, true, false, loadRead2, loadPairs);
    	if(!loadT1Reads || loadT2Reads)
			System.err.println("TophatFileHitLoader: You asked to load only Type1 or Type2 reads, we do not load this information from Tophat SAM format.");
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
		while (iter.hasNext()) {
		    SAMRecord record = iter.next();
		    
		    if (record.getReadUnmappedFlag()) {continue; }
		    if(record.getReadPairedFlag() && !record.getFirstOfPairFlag() && !loadRead2){continue;}
		    float weight = 1/(float)record.getIntegerAttribute("NH");
		    
		    Read currRead = new Read();
	        
		    List<AlignmentBlock> blocks = record.getAlignmentBlocks();
		    for(int a=0; a<blocks.size(); a++){ //Iterate over alignment blocks
		    	AlignmentBlock currBlock = blocks.get(a);
		    	int aStart = currBlock.getReferenceStart();
		    	int aEnd = aStart + currBlock.getLength()-1;
		    	int aLen = currBlock.getLength();
		    	boolean nearbyBlocks=true;
		    	while(nearbyBlocks && a<blocks.size()-1){
		    		if(blocks.get(a+1).getReferenceStart() - currBlock.getReferenceStart() < record.getReadLength()){
		    			aEnd = blocks.get(a+1).getReferenceStart() + blocks.get(a+1).getLength()-1;
		    			aLen += blocks.get(a+1).getLength();
		    			a++;
		    		}else{
		    			nearbyBlocks=false;
		    		}
		    	}
		    	
		    	ReadHit currHit = new ReadHit(
		    			record.getReferenceName().replaceFirst("^chr", ""),
		    			aStart, aEnd,
		    			record.getReadNegativeStrandFlag() ? '-' : '+',
		    			weight);
		   
		    	currRead.addHit(currHit);
			}	
		    addHits(currRead);
		    
		    //load pair if this is a first mate, congruent, proper pair
		    if(record.getFirstOfPairFlag() && record.getProperPairFlag()){
		    	boolean neg = record.getReadNegativeStrandFlag();
                boolean mateneg = record.getMateNegativeStrandFlag();
                HitPair hp = new HitPair((neg ? record.getAlignmentEnd() : record.getAlignmentStart()),
                		record.getMateReferenceName(),
                		(mateneg ? record.getMateAlignmentStart()+record.getReadLength()-1 : record.getMateAlignmentStart()), 
                		mateneg ? 1 : 0,
                		1);
                addPair(record.getReferenceName(), neg ? '-':'+', hp);
		    }
		}
		iter.close();
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }//end of sourceAllHits method
}
