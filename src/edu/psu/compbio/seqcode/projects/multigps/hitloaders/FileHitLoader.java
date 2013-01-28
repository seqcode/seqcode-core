package edu.psu.compbio.seqcode.projects.multigps.hitloaders;

import java.io.File;

/**
 * FileHitLoader: Loads reads from a collection of files. 
 * Formats supported:
 * ELAND, NOVO, BOWTIE, BED	, SAM, TOPSAM
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
	public FileHitLoader(File file, boolean useNonUnique){
		super();
		this.file = file;
		this.useNonUnique=useNonUnique;
	}
	
	/**
	 * No cleanup for file loaders
	 */
	public void cleanup(){}
}

