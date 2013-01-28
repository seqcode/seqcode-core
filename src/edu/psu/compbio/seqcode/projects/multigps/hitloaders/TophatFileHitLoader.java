package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.File;
import java.util.List;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import edu.psu.compbio.seqcode.projects.multigps.framework.*;

public class TophatFileHitLoader extends FileHitLoader{

    public TophatFileHitLoader(File f, boolean nonUnique) {
    	super(f, nonUnique);
    }

    /**
	 * Get the reads from the appropriate source (implementation-specific).
	 * Loads data to the fivePrimesList and hitsCountList
	 */
	public void sourceReads() {
		this.initialize();
		SAMFileReader reader = new SAMFileReader(file);
		reader.setValidationStringency(ValidationStringency.SILENT);
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
		reader.close();
    }//end of countReads method
}
