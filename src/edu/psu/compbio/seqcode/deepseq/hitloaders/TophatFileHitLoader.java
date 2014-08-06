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

import edu.psu.compbio.seqcode.deepseq.Read;
import edu.psu.compbio.seqcode.deepseq.ReadHit;

/**
 * TophatFileHitLoader: A FileHitLoader for Tophat SAM/BAM output.
 * Each alignment block is loaded as a separate hit. Ignores uniqueness. 
 * @author mahony
 *
 */
public class TophatFileHitLoader extends FileHitLoader{

    public TophatFileHitLoader(File f, boolean nonUnique) {
    	super(f, nonUnique);
    }

    /**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
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
		}
		iter.close();
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }//end of sourceAllHits method
}
