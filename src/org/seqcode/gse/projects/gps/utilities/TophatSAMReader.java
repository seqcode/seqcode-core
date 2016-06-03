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

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

public class TophatSAMReader extends AlignmentFileReader{

    public TophatSAMReader(File f, Genome g, int mis, boolean nonUnique, int idSeed) {
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
		while (iter.hasNext()) {
		    currID++;
		    SAMRecord record = iter.next();
		    
		    if (record.getReadUnmappedFlag()) {continue; }
		    float weight = 1/(float)record.getIntegerAttribute("NH");
		    if(readLength ==-1)
		    	readLength = record.getReadLength();
		    
		    Read currRead = new Read((int)totalWeight);
	        
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
		    	
		    	ReadHit currHit = new ReadHit(gen,
					  currID,
					  record.getReferenceName().replaceFirst("^chr", ""), 
					  aStart, aEnd, 
					  record.getReadNegativeStrandFlag() ? '-' : '+',
					  weight);
		   
		    	currRead.addHit(currHit);
		    	currID++;
			}	
		    addHits(currRead);
		    totalWeight++;
		}
		iter.close();
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		populateArrays();
    }//end of countReads method
    
}
