package org.seqcode.deepseq.hitloaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Arrays;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

import org.seqcode.genome.Genome;
import org.seqcode.deepseq.HitPair;
import org.seqcode.deepseq.Read;
import org.seqcode.deepseq.ReadHit;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.utils.HierarchicalHitInfo;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

public class HDF5HitLoader extends HitLoader {
	
	private boolean useChimericReads=false; //Ignore chimeric mappings for now.
	protected Genome genome;
	protected HierarchicalHitInfo readHHI;
	protected HierarchicalHitInfo pairHHI;
	
	// sam load identifier
	protected File file; 
	protected boolean isCache=false; //Read hits from a pre-stored HDF5HitCache hits info file
	protected boolean useNonUnique=false; //whether use non unique reads
	protected boolean loadRead2=true; //Load read 2 in paired-end
	protected boolean loadReads=true; //Load read information
	protected boolean loadPairs=false; //Load pair information (if exists)
	protected boolean hasPairs = false; //Flag to say there are pairs in the sample 
	protected int totalReads = 0; //total number of reads
	protected int totalHitPairs = 0; //total number of hit pairs
	protected String sourceName=""; //String describing the source
	
	protected LinkedBlockingQueue<ReadHit> readQueue = new LinkedBlockingQueue<ReadHit>(10000) ;
	protected LinkedBlockingQueue<HitPair> pairQueue = new LinkedBlockingQueue<HitPair>(10000) ;
	
	protected boolean terminatorFlag = false;
	
	public HDF5HitLoader(Genome g, File f, boolean loadReads, boolean nonUnique, boolean loadRead2, boolean loadPairs, boolean isCache) {
		super(true, false, loadRead2, loadPairs);
		this.genome = g;
		this.file = f;
		this.useNonUnique = nonUnique;
		this.loadReads = loadReads;
		this.loadPairs = loadPairs;
		this.sourceName = f.getAbsolutePath();
		this.isCache = isCache;
		
		if(!isCache) {
			this.readHHI = new HierarchicalHitInfo(genome, f.getAbsolutePath() + ".read.h5", false);
			this.pairHHI = new HierarchicalHitInfo(genome, f.getAbsolutePath() + ".pair.h5", true);
			readHHI.initializeHDF5();
			pairHHI.initializeHDF5();
		} else {
			this.readHHI = new HierarchicalHitInfo(genome, f.getAbsolutePath() + ".read.h5", false);
			this.pairHHI = new HierarchicalHitInfo(genome, f.getAbsolutePath() + ".pair.h5", true);
			readHHI.openDataSet();
			pairHHI.openDataSet();
		}
	}
	
	//Accessors
	public boolean isCache() {return isCache;}
	public boolean hasPairedReads() {return hasPairs;}
	public double getHitCount() {return totalReads;}
	public String getSourceName() {return sourceName;}
	public HierarchicalHitInfo getHitPairInfo() {return pairHHI;}
	public HierarchicalHitInfo getReadInfo() {return readHHI;}
	
	class ProcessPairThread implements Runnable {
		private ArrayList<HitPair> processingPairs;
		
		public void run() {
			while(!terminatorFlag || !pairQueue.isEmpty()) {
				processingPairs = new ArrayList<HitPair>();
				pairQueue.drainTo(processingPairs, 10000);
				process(processingPairs);
			}
		}
		
		private void process(List<HitPair> input) {
			for(HitPair hp : input) {
				if(!hasPairs) {
					hasPairs = true;
				}
				try {
					pairHHI.appendHit(hp.r1Chr, hp.r1Strand, new double[] {hp.r1Pos, genome.getChromID(hp.r2Chr), hp.r2Pos, hp.r2Strand, hp.pairWeight, hp.pairMid});
					totalHitPairs += (int)hp.pairWeight;
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
				readQueue.drainTo(processReads, 10000);
				process(processReads);
			}
		}
		
		private void process(List<ReadHit> input) {
			for(ReadHit rh : input) {
				try {
					readHHI.appendHit(rh.getChrom(), rh.getStrand()=='+' ? 0 : 1, new double[] {rh.getStrand()=='+' ? rh.getStart() : rh.getEnd(), rh.getStrand()=='+' ? rh.getEnd() : rh.getStart(), rh.getWeight()});
					totalReads += rh.getWeight();
				} catch (Exception e) {
					// TODO: handle exception
					e.printStackTrace();
				}
			}
		}
	}
		
	public void sourceAllHits() {
		if(isCache)
			return;
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
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
		    
		    if(record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary() && !useChimericReads){continue;}
		    if(record.getReadPairedFlag() && record.getSecondOfPairFlag() && !loadRead2){continue;}
		    
		    if(loadReads) {
			    if (lastread == null || !lastread.equals(record.getReadName())) {
			    	processRead(byRead);
			    	byRead.clear();
			    }
			    lastread = record.getReadName();
			    
			    byRead.add(record); //Filter by first or second of pair here if loading by type1/2?
		    }
		    
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
                    hp = null;
				} catch (Exception e) {
					// TODO: handle exception
					e.printStackTrace();
				}

		    }
		}
		terminatorFlag = true;
		try {
			ppt.join();
			prt.join();
			readQueue = null;
			pairQueue = null;
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
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
		
		currRead = null;
    }//end of processRead
	
    // close the hitloader
    public void close() {
    	readHHI.closeDataset();
    	readHHI.closeFile();
    	pairHHI.closeDataset();
    	pairHHI.closeFile();
    }
    
    // clean the dataset
    public void cleanup() {
    	if(!isCache) {
    		readHHI.deleteFile();
    		pairHHI.deleteFile();
    	}
    }
    
	public static void main(String[] args) {
		File genFile = new File(args[0]);
		File bamFile = new File(args[1]);
		
		Genome gen = new Genome("Genome", genFile, true);
		
		HDF5HitLoader hl = new HDF5HitLoader(gen, bamFile, false, false, false, true, false);
		
		hl.sourceAllHits();
		hl.getHitPairInfo().sortByReference();
		
		hl.close();
	}
	
}
