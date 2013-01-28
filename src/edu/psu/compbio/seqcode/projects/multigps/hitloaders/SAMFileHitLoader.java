package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import edu.psu.compbio.seqcode.projects.multigps.framework.*;

public class SAMFileHitLoader extends FileHitLoader{

    public SAMFileHitLoader(File f, boolean nonUnique) {
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
		Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
		String lastread = null;
		while (iter.hasNext()) {
		    SAMRecord record = iter.next();
		    
		    if (record.getReadUnmappedFlag()) {continue; }
		    if (lastread == null || !lastread.equals(record.getReadName())) {
		    	processRead(byRead);
		    	byRead.clear();
		    }
		    lastread = record.getReadName();
		    byRead.add(record);
			    
		}
		processRead(byRead);
		iter.close();
		reader.close();
    }//end of countReads method
    
    protected void processRead(Collection<SAMRecord> records) {
        int mapcount = records.size();
        if (mapcount == 0) {
            return;
        }
        if (!useNonUnique && mapcount > 1) {
            return;
        }
	
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
