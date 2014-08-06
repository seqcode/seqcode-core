package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.SeqBiasModel;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;


public class SAMReadLoader extends ReadLoader{

	protected File samFile;
	protected SamReader reader=null;
	protected CloseableIterator<SAMRecord> readIterator=null;
	protected SAMRecord currRecord;
	
	public SAMReadLoader(File samFile, String cond, Genome gen, Map<String, String> nameTrans){
		super(gen, nameTrans);
		this.samFile = samFile;
		sourcePath = samFile.getAbsolutePath();
		sourceName = samFile.getName();
		conditionName = cond;
		SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		reader = factory.open(samFile);
	}

	
	public Genome estimateGenome() {
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();
		if(dictionary !=null){
			for(SAMSequenceRecord record : dictionary.getSequences()){
			    String chr = record.getSequenceName();
			    chrLenMap.put(chr, record.getSequenceLength());
			}
		}else{
			CloseableIterator<SAMRecord> iter = reader.iterator();
			while (iter.hasNext()) {
			    SAMRecord record = iter.next();
			    String chr = record.getReferenceName();
			    int max = Math.max(record.getAlignmentEnd(), record.getAlignmentStart());
			    if(!chrLenMap.containsKey(chr) || chrLenMap.get(chr)<max)
					chrLenMap.put(chr, max);
			}iter.close();
		}
		return(new Genome("Genome", chrLenMap));
	}

	/**
	 * Get reads that overlap a region. Uses picard search method, so can only be used for BAM files. 
	 * All other SAM Iterators need to be closed before calling this method.
	 * @param r Region to query
	 * @return A List of AlignHits that overlap the Region 
	 */
	public List<AlignHit> getOverlappingHits(Region r, SeqBiasModel bias) {
		List<AlignHit> hits = new ArrayList<AlignHit>();
		
		SAMRecordIterator iter = reader.queryOverlapping(chromNameTranslator.get(r.getChrom()), r.getStart(), r.getEnd());
		while(iter.hasNext()){
			SAMRecord s  = iter.next();
			double w = bias==null ? 1.0 : bias.getWeight(s);
			hits.add(new AlignHit(s, gen, w));
		}
		iter.close();
		
		return hits;
	}
	
	/**
	 * Get reads that overlap a region from the stream. 
	 * Assumes that reads are sorted, and that regions are queried in the same order. 
	 * This method should only be used to get reads for large non-overlapping regions. If overlapping regions are queried (or neighboring exons), you will miss out on some reads.
	 * @param r Region to query
	 * @return A List of AlignHits that overlap the Region 
	 */
	public List<AlignHit> getOverlappingHitsFromStream(Region r, SeqBiasModel bias) {
		List<AlignHit> hits = new ArrayList<AlignHit>();
		
		while(currRecord.getReferenceName().compareTo(chromNameTranslator.get(r.getChrom()))<0 && nextRead()){}
		
		boolean moreReads = readIterator.hasNext(); 
		while(currRecord.getReferenceName().compareTo(chromNameTranslator.get(r.getChrom()))==0 && currRecord.getAlignmentStart()<r.getEnd() && moreReads){
			//The overlap comparison here assumes that we are on the same chromosome and that the alignment start can't be greater than the region end
			if(currRecord.getAlignmentStart()>r.getStart() || (currRecord.getAlignmentEnd()>r.getStart()&&currRecord.getAlignmentEnd()<r.getEnd())){
				double w = bias==null ? 1.0 : bias.getWeight(currRecord);
				hits.add(new AlignHit(currRecord, gen, w));
			}
			moreReads = nextRead();
		}
		
		return hits;
	}

	public void initializeReadIterator() {
		closeReader();
		//Reopen SAM file. Seems to have been closed when closing the read iterator
		SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		reader = factory.open(samFile);
		readIterator = reader.iterator(); 
		nextRead();
	}
	
	public void closeReader() {
		if(readIterator!=null){
			readIterator.close();
			readIterator=null;
		}
		if(reader != null){
			try {
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public boolean iteratorHasNextHit(){
		if(readIterator.hasNext())
			return true;
		else
			return false;
	}
	
	public boolean nextRead() {
		if(readIterator.hasNext()){
			currRecord = readIterator.next();
			return true;
		}else
			return false;
	}

	public AlignHit getCurrHit() {
		return (new AlignHit(currRecord, gen));
	}

	public byte[] getCurrReadSeq() {
		return(currRecord.getReadBases());
	}
}

