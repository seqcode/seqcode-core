package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.SeqBiasModel;

public class SAMReadLoader extends ReadLoader{

	protected File samFile;
	protected SAMFileReader reader=null;
	protected CloseableIterator<SAMRecord> readIterator=null;
	protected SAMRecord currRecord;
	
	public SAMReadLoader(File samFile, String cond, Genome gen, Map<String, String> nameTrans){
		super(gen, nameTrans);
		this.samFile = samFile;
		sourcePath = samFile.getAbsolutePath();
		sourceName = samFile.getName();
		conditionName = cond;
		reader = new SAMFileReader(samFile);
		reader.setValidationStringency(ValidationStringency.SILENT);
		reader.enableIndexMemoryMapping(false);
		reader.enableIndexCaching(false);
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
		reader = new SAMFileReader(samFile);
		reader.setValidationStringency(ValidationStringency.SILENT);
		readIterator = reader.iterator(); 
		nextRead();
	}
	
	public void closeReader() {
		if(readIterator!=null){
			readIterator.close();
			readIterator=null;
		}
		if(reader != null){
			reader.close();
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

