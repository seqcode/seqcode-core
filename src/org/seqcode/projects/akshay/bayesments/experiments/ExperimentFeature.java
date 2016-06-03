package org.seqcode.projects.akshay.bayesments.experiments;

import java.util.*;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.projects.akshay.bayesments.framework.BayesmentsConfig;


public class ExperimentFeature {
	
	protected int index;
	protected BayesmentsConfig config;
	protected List<Sample> signalSamples = new ArrayList<Sample>();
	protected List<Sample> controlSamples = new ArrayList<Sample>();
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected List<ExperimentCondition> conditions = new ArrayList<ExperimentCondition>();
	protected HashMap<Integer,ExperimentCondition> indexedConditions = new HashMap<Integer, ExperimentCondition>();
	protected HashMap<ExperimentCondition, Integer> conditionIndex = new HashMap<ExperimentCondition, Integer>();
	protected String name;
	
	
	public ExperimentFeature(BayesmentsConfig config, int id, String name, List<ExperimentCondition> conditions) {
		this.index = id;
		this.config = config;
		this.name = name;
		this.conditions = conditions;
		
		for(ExperimentCondition ec : conditions){
			List<ControlledExperiment> tempRepList = ec.getReplicates();
			for (ControlledExperiment ce : tempRepList){
				if(!replicates.contains(ce)){
					this.replicates.add(ce);
				}
			}
			this.indexedConditions.put(ec.getIndex(), ec);
			this.conditionIndex.put(ec, ec.getIndex());
		}
		
		
		
		for(ControlledExperiment ce : replicates){
			if(!signalSamples.contains(ce.getSignal())){
				signalSamples.add(ce.getSignal());
			}
			if(!controlSamples.contains(ce.getControl())){
				controlSamples.add(ce.getControl());
			}
		}
	}
	
	
	//Accessors
	
	public int getIndex(){return this.index;}
	public String getName(){return this.name;}
	public List<ExperimentCondition> getCondtionList(){return this.conditions;}
	public int getNumReplicates(){return this.replicates.size();}
	public int getNumCondition(){return this.conditions.size();}
	public String getNameOfConditionAtGivenIndexInConditionList(int i){return conditions.get(i).getName();}
	
}
