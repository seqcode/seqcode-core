package edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments;

import java.util.*;

import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;

public class ExperimentSet {
	
	protected List<Sample> samples = new ArrayList<Sample>();
	protected List<ExperimentCondition> conditions =  new ArrayList<ExperimentCondition>();
	protected List<ControlledExperiment> replicates =  new ArrayList<ControlledExperiment>();
	protected List<ExperimentFeature> features = new ArrayList<ExperimentFeature>();
	protected HashMap<String, ExperimentCondition> namedConditions = new HashMap<String, ExperimentCondition>();
	protected HashMap<Integer, ExperimentCondition> indexedConditions =  new HashMap<Integer, ExperimentCondition>();
	protected HashMap<ExperimentCondition, Integer> conditionIndex =  new HashMap<ExperimentCondition, Integer>();
	protected HashMap<String, ExperimentFeature> namedFeatures = new HashMap<String, ExperimentFeature>();
	protected HashMap<Integer, ExperimentFeature> indexedFeatures = new HashMap<Integer, ExperimentFeature>();
	protected HashMap<ExperimentFeature, Integer> featureIndex = new HashMap<ExperimentFeature, Integer>();
	
	public ExperimentSet(List<ControlledExperiment> replicates, List<ExperimentCondition> conditions, List<ExperimentFeature> features) {
		this.replicates = replicates;
		this.conditions = conditions;
		this.features = features;
		
		for(ControlledExperiment ce : this.replicates){
			if(!samples.contains(ce.getSignal())){
				this.samples.add(ce.getSignal());
			}
			if(!samples.contains(ce.getControl())){
				this.samples.add(ce.getControl());
			}
		}
		
		for(ExperimentCondition con : this.conditions){
			this.namedConditions.put(con.getName(), con);
			this.conditionIndex.put(con, con.getIndex());
			this.indexedConditions.put(con.getIndex(), con);
		}
		
		for(ExperimentFeature fea : this.features){
			this.namedFeatures.put(fea.getName(), fea);
			this.featureIndex.put(fea, fea.getIndex());
			this.indexedFeatures.put(fea.getIndex(), fea);
		}
		
	}
	
	//Accessors
	
	public int getNumfChromatinFeatures(){return this.namedFeatures.get("CHROMATIN").getNumReplicates();}
	
	

}
