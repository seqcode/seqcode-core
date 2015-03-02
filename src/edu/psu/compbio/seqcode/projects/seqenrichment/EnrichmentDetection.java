package edu.psu.compbio.seqcode.projects.seqenrichment;

import java.util.List;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;


/**
 * EnrichmentDetection is the parent class for all enrichment-based peak/domain callers in this package. 
 *  
 * Handles genome & experiment & configuration loading. 
 * 
 * @author shaunmahony
 *
 */
public abstract class EnrichmentDetection {
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	
	//Constructors
	public EnrichmentDetection(String []args){
		gconfig = new GenomeConfig(args);
		econfig = new ExptConfig(gconfig.getGenome(), args);
		manager = new ExperimentManager(econfig);
		
	}
	
	//public abstract List<Feature> execute();	
	public abstract void printError();
}
