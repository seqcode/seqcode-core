package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.util.ArrayList;

/**
 * ExptCondition: an object representing a collection of experiment replicates from the same experimental condition
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExptCondition {

	protected String name;
	protected ArrayList<ExptReplicate> replicates;
	
	public ExptCondition(String name){
		this.name = name;
		replicates = new ArrayList<ExptReplicate>();
	}
	
	//Accessors
	public String getName(){return name;}
	public ArrayList<ExptReplicate> getReplicates(){return replicates;}
	
	/**
	 * Add a new replicate to the condition
	 * @param rep Instantiated ExptReplicate
	 */
	public void addReplicate(ExptReplicate rep){
		replicates.add(rep);
	}
	
	/**
	 * Reset the read loaders
	 */
	public void resetHitExtractors(){
		for(ExptReplicate r : replicates){
			r.resetHitExtractor();
		}
	}
}
