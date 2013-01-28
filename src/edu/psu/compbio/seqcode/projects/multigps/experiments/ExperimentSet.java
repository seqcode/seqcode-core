package edu.psu.compbio.seqcode.projects.multigps.experiments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


/**
 * ExperimentSet is a collection of ExperimentConditions. 
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentSet {

	protected List<Sample> samples = new ArrayList<Sample>();
	protected List<ExperimentCondition> conditions = new ArrayList<ExperimentCondition>();
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected HashMap<ExperimentCondition, Integer> conditionIndex = new HashMap<ExperimentCondition, Integer>();
	protected HashMap<Integer, ExperimentCondition> indexedCondition = new HashMap<Integer, ExperimentCondition>();
	protected HashMap<String, ExperimentCondition> namedCondition = new HashMap<String, ExperimentCondition>();
	protected int numConditions = 0;
	
	public ExperimentSet(List<ExperimentCondition> conds, List<ControlledExperiment> reps){
		conditions = conds;
		replicates = reps;
		numConditions = conditions.size();
		
		//Fill Sample list
		for(ControlledExperiment c : replicates){
			if(!samples.contains(c.getSignal()))
				samples.add(c.getSignal());
			if(!samples.contains(c.getControl()))
				samples.add(c.getControl());
		}
		
		for(int i=0; i<numConditions; i++){
			conditionIndex.put(conditions.get(i), i);
			indexedCondition.put(i, conditions.get(i));
			namedCondition.put(conditions.get(i).getName(), conditions.get(i));
			//System.err.println("Condition index check: "+conditions.get(i).getName()+" "+conditions.get(i).getIndex()+" "+i);
		}
	}
	
	//Accessors
	public List<Sample> getSamples(){return samples;}
	public List<ExperimentCondition> getConditions(){return conditions;}
	public List<ControlledExperiment> getReplicates(){return replicates;}
	public int getConditionIndex(ExperimentCondition c){return conditionIndex.get(c);}
	public ExperimentCondition getIndexedCondition(int index){return indexedCondition.get(index);}
	public ExperimentCondition getNamedCondition(String name){return namedCondition.get(name);}
}
