package edu.psu.compbio.seqcode.projects.seqenrichment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.projects.seqenrichment.features.Feature;


/**
 * FeatureDetection is the parent class for all enrichment-based peak/domain callers in this package. 
 *  
 * Handles genome & experiment & configuration loading. 
 * 
 * @author shaunmahony
 *
 */
public abstract class FeatureDetection {
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	protected SeqEnrichmentConfig seconfig;
	protected Map<ExperimentCondition, List<Feature>> features;
	
	//Constructors
	public FeatureDetection(String []args){
		gconfig = new GenomeConfig(args);
		econfig = new ExptConfig(gconfig.getGenome(), args);
		seconfig = new SeqEnrichmentConfig(args);
		manager = new ExperimentManager(econfig);
		
		features=new HashMap<ExperimentCondition, List<Feature>>();
	}
	
	public abstract Map<ExperimentCondition, List<Feature>> execute();
	public abstract void printError();
	
	public List<Feature> getFeatures(ExperimentCondition c){
		if(features==null)
			return null;
		else if(!features.containsKey(c))
			return new ArrayList();
		else
			return features.get(c);
	}
}
