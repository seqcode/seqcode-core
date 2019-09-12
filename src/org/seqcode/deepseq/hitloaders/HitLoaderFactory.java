package org.seqcode.deepseq.hitloaders;

import java.io.File;
import java.io.ObjectInputFilter.Config;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.seqcode.data.seqdata.SeqDataLoader;
import org.seqcode.data.seqdata.SeqLocator;
import org.seqcode.deepseq.experiments.ExptConfig;

/**
 * HitLoaderFactory is a simple class that defines the hitloaders. 
 * @author mahony
 *
 */
public class HitLoaderFactory {

	ExptConfig econfig; 
	
	public HitLoaderFactory(ExptConfig e){
		econfig = e;
	}
	
	/**
	 * Add a ReadDB HitLoader.
	 * 
	 * @param locs List of SeqLocators
	 */
	public HitLoader makeReadDBHitLoader(SeqDataLoader loader, String name){
		if(loader==null)
			throw new RuntimeException("SeqDataLoader not initialized in makeReadDBHitLoader");
		
		List<SeqLocator> locs = new ArrayList<SeqLocator>();
		String[] pieces = name.trim().split(";");
		Set<String> reps = new TreeSet<String>(); 
		if (pieces.length == 2) {
			//Merge all available replicates
			locs.add(new SeqLocator(pieces[0], reps, pieces[1]));
        } else if (pieces.length == 3) {
        	//Specific named replicate or merge a list of comma-separated replicates
        	String[] rs = pieces[1].split(",");  
        	for(int r=0; r<rs.length; r++)
        		reps.add(rs[r]);
        	locs.add(new SeqLocator(pieces[0], reps, pieces[2]));
        } else {
            throw new RuntimeException("Couldn't parse a SeqLocator from " + name);
        }
		return (new ReadDBHitLoader(loader, econfig.getGenome(), locs, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs()));
	}
	

	/**
	 * Add a File HitLoader. File formats accepted include:
	 * SCIDX, NOVO, BOWTIE, BED, SAM, TOPSAM
	 * @param files List of File/String Pairs, where the string is a format descriptor
	 */
	public HitLoader makeFileHitLoader(String filename, String format, boolean useNonUnique){
		HitLoader currReader=null;
		File file = new File(filename);
		if(!file.isFile() && !econfig.getKeepHDF5()){System.err.println("File not found: "+file.getName());System.exit(1);}
		if(format.equals("SAM") || format.equals("BAM")){
			currReader = new SAMFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadRead2(), econfig.getLoadPairs());
		}else if(format.equals("TOPSAM")){
			currReader = new TophatFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(),  econfig.getLoadRead2(), econfig.getLoadPairs());
		}else if(format.equals("NOVO")){
			currReader = new NovoFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("BOWTIE")){
			currReader = new BowtieFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("BED")){
			currReader = new BEDFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("SCIDX") || format.equals("IDX")){
			currReader = new IDXFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("HDF5")) {
			currReader = new HDF5HitLoader(econfig.getGenome(), file, econfig.getLoadReads(), useNonUnique, econfig.getLoadRead2(), econfig.getLoadPairs(), false);			
		}else if(format.equals("HDF5Cache")) {
			currReader = new HDF5HitLoader(econfig.getGenome(), file, econfig.getLoadReads(), useNonUnique, econfig.getLoadRead2(), econfig.getLoadPairs(), true);			
		}
		else{
		    System.err.println("Unknown file format: "+format);
		    System.exit(1);
		}
		return currReader;
	}
}
