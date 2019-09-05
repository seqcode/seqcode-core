package org.seqcode.deepseq.hitloaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import org.seqcode.genome.Genome;
import org.seqcode.deepseq.HitPair;

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
	
	// sam load identifier
	protected File file; 
	protected boolean loadType1=true; //Load type1 reads
	protected boolean loadType2=false; //Load type2 reads (if exists)
	protected boolean loadRead2=true; //Load read 2 in paired-end
	protected boolean loadPairs=false; //Load pair information (if exists)
	protected boolean hasPairs = false; //Flag to say there are pairs in the sample 
	protected double totalHits; //totalHits is the sum of alignment weights
	protected String sourceName=""; //String describing the source
	
	// hdf5 identifier
	protected long file_id = -1;
	protected long group_id = -1;
	protected long dataspace_id = -1;
	protected long dataset_id = -1;
	protected long dcpl_id = -1;
	protected String GROUPNAME = "group";
	protected String DATASETNAME = "dataset";
	
	protected int nChroms = 15;
	protected LinkedBlockingQueue<HitPair> readQueue = new LinkedBlockingQueue<HitPair>() ;
	
	public HDF5HitLoader(Genome g, File f, boolean nonUnique, boolean loadT1Reads, boolean loadT2Reads, boolean loadRead2, boolean loadPairs) {
		this.genome = g;
		this.file = f;
		this.loadType1 = loadT1Reads;
		this.loadType2 = loadT2Reads;
		this.loadRead2 = loadRead2;
		this.loadPairs = loadPairs;
		totalHits = 0;
		
	}
	
	class ProcessReadThread implements Runnable {
		private boolean terminatorFlag = false;
		private List<HitPair> processingReads = new ArrayList<HitPair>();
		
		public void run() {
			while(!terminatorFlag) {
				readQueue.drainTo(processingReads, 100);
				process(processingReads);
			}
		}
		
		private void process(List<HitPair> input) {
			for(HitPair hp : input) {
				int dim1 = convertIndex(hp.r1Chr, hp.r1Strand);
			}
		}
		
		private int convertIndex(String chr, int strand) {
			return genome.getChromID(chr) + strand;
		}
	}
	
	public void initializeHDF5() {
		// Create a new file using default properties.
		try {
			file_id = H5.H5Fcreate("test.hdf5", HDF5Constants.H5F_ACC_TRUNC, 
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Create a new group
		try {
			if(file_id >= 0)
				group_id = H5.H5Gcreate(file_id, "/" + GROUPNAME, HDF5Constants.H5P_DEFAULT,
						 HDF5Constants.H5P_DEFAULT,  HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Create the dataspace
		try {
			dataspace_id = H5.H5Screate_simple(2, new long[] {nChroms, 10000}, null);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Create the dataset creation property list
		try {
			dcpl_id = H5.H5Pcreate(HDF5Constants.H5P_DATASET_CREATE);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Set the chunk size
		try {
			if (dcpl_id >= 0) 
				H5.H5Pset_chunk(dcpl_id, 1, new long[] {nChroms, 10});
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Create the chunked dataset
		try {
			dataset_id = H5.H5Dcreate(file_id, "/" + GROUPNAME + "/" + DATASETNAME, HDF5Constants.H5T_STD_I32BE, 
					dataspace_id, HDF5Constants.H5P_DEFAULT, dcpl_id,HDF5Constants.H5P_DEFAULT);
		} catch (Exception e) {
			e.printStackTrace();
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
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
		    
		    if(record.getReadUnmappedFlag()) {continue; }
		    if(record.isSecondaryOrSupplementary() && !useChimericReads){continue;}
		    if(record.getReadPairedFlag() && record.getSecondOfPairFlag() && !loadRead2){continue;}
		    	
		    //load pair if this is a first mate, congruent, proper pair
		    if(loadPairs && record.getFirstOfPairFlag() && record.getProperPairFlag()){
		    	boolean neg = record.getReadNegativeStrandFlag();
                boolean mateneg = record.getMateNegativeStrandFlag();
                HitPair hp = new HitPair((neg ? record.getAlignmentEnd() : record.getAlignmentStart()),
                		record.getMateReferenceName().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr", ""),
                		(mateneg ? record.getMateAlignmentStart()+record.getReadLength()-1 : record.getMateAlignmentStart()), 
                		mateneg ? 1 : 0,
                		1);
                readQueue.add(hp);
		    }
		}
	}
}
