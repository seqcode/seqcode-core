package org.seqcode.gse.projects.gps.utilities;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.gse.projects.gps.Read;
import org.seqcode.gse.projects.gps.ReadHit;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

public class SAMReader extends AlignmentFileReader{

    public SAMReader(File f, Genome g, int mis, boolean nonUnique, int idSeed) {
    	super(f, g, mis, nonUnique, idSeed);
    }
    
	protected void estimateGenome() {
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		SamReader reader = factory.open(inFile);
		SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();
		if(dictionary !=null){
			for(SAMSequenceRecord record : dictionary.getSequences()){
			    String chr = record.getSequenceName().replaceFirst("^chr", "");
			    chrLenMap.put(chr, record.getSequenceLength());
			}
		}else{
			CloseableIterator<SAMRecord> iter = reader.iterator();
			while (iter.hasNext()) {
			    currID++;
			    SAMRecord record = iter.next();
			    String chr = record.getReferenceName().replaceFirst("^chr", "");
			    int max = Math.max(record.getAlignmentEnd(), record.getAlignmentStart());
			    if(!chrLenMap.containsKey(chr) || chrLenMap.get(chr)<max)
					chrLenMap.put(chr, max);
			}
		}
		gen=new Genome("Genome", chrLenMap);
	}
	
    //Return the total reads and weight
    protected void countReads() {
		readLength=-1;
		totalHits=0;
		totalWeight=0;
		
		SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		SamReader reader = factory.open(inFile);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
		String lastread = null;
		while (iter.hasNext()) {
		    currID++;
		    SAMRecord record = iter.next();
		    if(readLength ==-1)
		    	readLength = record.getReadLength();
		    
		    if (record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary()){continue;}
		    if (lastread == null || !lastread.equals(record.getReadName())) {
		    	processRead(byRead);
		    	byRead.clear();
		    }
		    lastread = record.getReadName();
		    byRead.add(record);
			    
		}
		processRead(byRead);
		iter.close();
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		populateArrays();
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
        Read currRead = new Read((int)totalWeight);
		for (SAMRecord record : records) {
			int start = record.getAlignmentStart();
		    int end = record.getAlignmentEnd();
		    ReadHit currHit = new ReadHit(gen, currID,
						  record.getReferenceName().replaceFirst("^chr", ""), 
						  start, end, 
						  record.getReadNegativeStrandFlag() ? '-' : '+',
						  weight);
		    currRead.addHit(currHit);
		    currID++;
		}currRead.setNumHits(mapcount);
		addHits(currRead);
		totalWeight++;
    }//end of processRead

}
