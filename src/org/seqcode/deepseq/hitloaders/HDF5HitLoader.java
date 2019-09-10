package org.seqcode.deepseq.hitloaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

import org.seqcode.genome.Genome;
import org.seqcode.deepseq.HitPair;
import org.seqcode.deepseq.Read;
import org.seqcode.deepseq.ReadHit;
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
	protected HierarchicalHitInfo readHHI;
	protected HierarchicalHitInfo pairHHI;
	
	// sam load identifier
	protected File file; 
	protected boolean useNonUnique=false; //whether use non unique reads
	protected boolean loadType1=true; //Load type1 reads
	protected boolean loadType2=false; //Load type2 reads (if exists)
	protected boolean loadRead2=true; //Load read 2 in paired-end
	protected boolean loadPairs=false; //Load pair information (if exists)
	protected boolean hasPairs = false; //Flag to say there are pairs in the sample 
	protected int totalReads = 0; //total number of reads
	protected int totalHitPairs = 0; //total number of hit pairs
	protected String sourceName=""; //String describing the source
	
	protected LinkedBlockingQueue<ReadHit> readQueue = new LinkedBlockingQueue<ReadHit>(1000) ;
	protected LinkedBlockingQueue<HitPair> pairQueue = new LinkedBlockingQueue<HitPair>(1000) ;
	
	protected boolean terminatorFlag = false;
	
	public HDF5HitLoader(Genome g, File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadRead2, boolean loadPairs) {
		this.genome = g;
		this.file = f;
		this.useNonUnique = nonUnique;
		this.loadType1 = loadT1Reads;
		this.loadType2 = loadT2Reads;
		this.loadRead2 = loadRead2;
		this.loadPairs = loadPairs;
		
		this.readHHI = new HierarchicalHitInfo(genome, f.getAbsolutePath() + ".read", false);
		this.pairHHI = new HierarchicalHitInfo(genome, f.getAbsolutePath() + ".pair", true);
		readHHI.initializeHDF5();
		pairHHI.initializeHDF5();
	}
	
	class ProcessPairThread implements Runnable {
		private ArrayList<HitPair> processingPairs;
		
		public void run() {
			while(!terminatorFlag || !pairQueue.isEmpty()) {
				processingPairs = new ArrayList<HitPair>();
				pairQueue.drainTo(processingPairs, 1000000);
				process(processingPairs);
			}
		}
		
		private void process(List<HitPair> input) {
			for(HitPair hp : input) {
				try {
					pairHHI.appendHit(hp.r1Chr, hp.r1Strand, new double[] {hp.r1Pos, genome.getChromID(hp.r2Chr), hp.r2Strand, hp.r2Pos, hp.pairWeight, hp.pairMid});
					totalHitPairs += 1;
				} catch (Exception e) {
					// TODO: handle exception
					e.printStackTrace();
				}
			}
		}
	}
	
	class ProcessReadThread implements Runnable {
		private ArrayList<ReadHit> processReads;
			
		public void run() {
			while (!terminatorFlag || !readQueue.isEmpty()) {
				processReads = new ArrayList<ReadHit>();
				readQueue.drainTo(processReads, 1000000);
				process(processReads);
			}
		}
		
		private void process(List<ReadHit> input) {
			for(ReadHit rh : input) {
				try {
					readHHI.appendHit(rh.getChrom(), rh.getStrand()=='+' ? 0 : 1, new double[] {rh.getStart(), rh.getEnd(), rh.getWeight()});
					totalReads += 1;
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
		Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
		String lastread = null;	
		// start threads used to handle input reads and hitpairs
		Thread ppt = new Thread(new ProcessPairThread());
		Thread prt = new Thread(new ProcessReadThread());
		prt.start();
		ppt.start();
		int count = 0;
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
		    
		    if(record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary() && !useChimericReads){continue;}
		    if(record.getReadPairedFlag() && record.getSecondOfPairFlag() && !loadRead2){continue;}
		    
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
                HitPair hp = new HitPair(record.getReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""),
                		(neg ? record.getAlignmentEnd() : record.getAlignmentStart()),
                		neg ? 1: 0,
                		record.getMateReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""),
                		(mateneg ? record.getMateAlignmentStart()+record.getReadLength()-1 : record.getMateAlignmentStart()), 
                		mateneg ? 1 : 0,
                		1);
                try {
                    pairQueue.put(hp);
				} catch (Exception e) {
					// TODO: handle exception
					e.printStackTrace();
				}

		    }
		    count++;
		    if(count % 1000000 == 0)
		    	System.err.format("Having loaded read: %d\t hitpair: %d\n", totalReads, totalHitPairs);
		}
		System.err.println("All reads have been loaded! loaded reads: " + totalReads + "\thitpairs: " + totalHitPairs);
		terminatorFlag = true;
		try {
			ppt.join();
			prt.join();
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		
		// print the length of each dimension
		for (String chr: genome.getChromList()) {
			System.out.println(pairHHI.getLength(chr, 0));
			System.out.println(pairHHI.getLength(chr, 1));
			System.out.println(readHHI.getLength(chr, 0));
			System.out.println(readHHI.getLength(chr, 1));
		}

		// close the dataset
		readHHI.closeDataset();
		readHHI.closeFile();
		pairHHI.closeDataset();
		pairHHI.closeFile();
	}
	
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
						  record.getReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""), 
						  start, end, 
						  record.getReadNegativeStrandFlag() ? '-' : '+',
						  weight);
		    currRead.addHit(currHit);
		}currRead.setNumHits(mapcount);

		for (ReadHit rh : currRead.getHits()) {
			try {
				readQueue.put(rh);
			} catch (Exception e) {
				// TODO: handle exception
				e.printStackTrace();
			}
		}
    }//end of processRead
	
	public static void main(String[] args) {
		File genFile = new File(args[0]);
		File bamFile = new File(args[1]);
		
		Genome gen = new Genome("Genome", genFile, true);
		
		for(String chr: gen.getChromList()) {
			System.out.println(chr + "\t" + gen.getChromID(chr));
		}
		HDF5HitLoader hl = new HDF5HitLoader(gen, bamFile, false, true, false, true, true);
		
		hl.sourceAllHits();
		
	}
	
}
