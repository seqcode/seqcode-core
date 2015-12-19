package edu.psu.compbio.seqcode.deepseq.hitloaders;

import java.io.File;

/**
 * FileHitLoader: Loads reads from a collection of files. 
 * Formats supported:
 * BOWTIE, BED, SAM, TOPSAM, NOVO
 * 
 * @author shaun
 *
 */
public abstract class FileHitLoader extends HitLoader{

	protected File file;
	protected boolean useNonUnique=true;
		
	/**
	 * Constructor
	 * @param g Genome
	 * @param name String
	 * @param files Pairs of Files and Strings (formats)
	 * @param useNonUnique boolean -- load non-uniquely mapping reads
	 */
	public FileHitLoader(File file, boolean useNonUnique, boolean loadT1, boolean loadT2, boolean loadRead2, boolean loadPairs){
		super(loadT1, loadT2, loadRead2, loadPairs);
		this.file = file;
		this.useNonUnique=useNonUnique;
		this.sourceName = file.getName();
	}
	
	/**
	 * No cleanup for file loaders
	 */
	public void cleanup(){}
}

