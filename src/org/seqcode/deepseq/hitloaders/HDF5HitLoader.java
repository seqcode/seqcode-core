package org.seqcode.deepseq.hitloaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import org.seqcode.genome.Genome;
import org.seqcode.deepseq.HitPair;
import org.seqcode.deepseq.utils.HierarchicalHitInfo;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

public class HDF5HitLoader {
	
	private boolean useChimericReads=false; //Ignore chimeric mappings for now.
	protected Genome genome;
	protected HierarchicalHitInfo hhInfo;
	
	// sam load identifier
	protected File file; 
	protected boolean loadType1=true; //Load type1 reads
	protected boolean loadType2=false; //Load type2 reads (if exists)
	protected boolean loadRead2=true; //Load read 2 in paired-end
	protected boolean loadPairs=false; //Load pair information (if exists)
	protected boolean hasPairs = false; //Flag to say there are pairs in the sample 
	protected double totalHits; //totalHits is the sum of alignment weights
	protected String sourceName=""; //String describing the source
	
	protected LinkedBlockingQueue<HitPair> readQueue = new LinkedBlockingQueue<HitPair>() ;
	
	protected boolean terminatorFlag = false;
	
	public HDF5HitLoader(Genome g, File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadRead2, boolean loadPairs) {
		this.genome = g;
		this.file = f;
		this.loadType1 = loadT1Reads;
		this.loadType2 = loadT2Reads;
		this.loadRead2 = loadRead2;
		this.loadPairs = loadPairs;
		totalHits = 0;
		
		this.hhInfo = new HierarchicalHitInfo(genome, f.getName(), 7);
		hhInfo.initializeHDF5();
	}
	
	class ProcessReadThread implements Runnable {
		private List<HitPair> processingReads = new ArrayList<HitPair>();
		
		public void run() {
			while(!terminatorFlag || !readQueue.isEmpty()) {
				readQueue.drainTo(processingReads, 100);
				process(processingReads);
			}
		}
		
		private void process(List<HitPair> input) {
			for(HitPair hp : input) {
				try {
					hhInfo.appendHit(hp.r1Chr, hp.r1Strand, new double[] {hp.r1Pos, genome.getChromID(hp.r2Chr), hp.r2Strand, hp.r2Pos, hp.pairWeight});
					totalHits += hp.pairWeight;
				} catch (Exception e) {
					// TODO: handle exception
					e.printStackTrace();
				}
			}
		}
	}
		
	public void sourceAllHits() {
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
				.validationStringency(ValidationStringency.SILENT);
		SamReader reader = factory.open(file);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		Thread prt = new Thread(new ProcessReadThread());
		prt.start();
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
		    
		    if(record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary() && !useChimericReads){continue;}
		    if(record.getReadPairedFlag() && record.getSecondOfPairFlag() && !loadRead2){continue;}
		    	
		    //load pair if this is a first mate, congruent, proper pair
		    if(loadPairs && record.getFirstOfPairFlag() && record.getProperPairFlag()){
		    	boolean neg = record.getReadNegativeStrandFlag();
                boolean mateneg = record.getMateNegativeStrandFlag();
                HitPair hp = new HitPair(record.getReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""),
                		neg ? 1: 0,
                		(neg ? record.getAlignmentEnd() : record.getAlignmentStart()),
                		record.getMateReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""),
                		(mateneg ? record.getMateAlignmentStart()+record.getReadLength()-1 : record.getMateAlignmentStart()), 
                		mateneg ? 1 : 0,
                		1);
                readQueue.add(hp);
		    }
		}
		terminatorFlag = true;
	}
	
	public static void main(String[] args) {
		File genFile = new File("D:\\Dropbox\\Code\\sem-test\\yeast\\geninfo.txt");
		File bamFile = new File("D:\\Dropbox\\Code\\sem-test\\yeast\\ref.F1804.chrI.bam");
		
		Genome gen = new Genome("Genome", genFile, true);
		HDF5HitLoader hl = new HDF5HitLoader(gen, bamFile, false, true, false, true, true);
		
		hl.sourceAllHits();
		
	}
	
}
