package edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments;

import java.util.ArrayList;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;

public class ExperimentFeature {
	
	protected int index;
	protected Config config;
	protected List<Sample> signalSamples = new ArrayList<Sample>();
	protected List<Sample> controlSamples = new ArrayList<Sample>();
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected HashMap<Integer,ControlledExperiment> indexedReplicates = new HashMap<Integer, ControlledExperiment>();
	protected HashMap<ControlledExperiment, Integer> replicateIndex = new HashMap<ControlledExperiment, Integer>();
	protected String name;
	
	
	public ExperimentFeature(Config config, int id, String name, List<ControlledExperiment> replicates) {
		this.index = id;
		this.config = config;
		this.name = name;
		this.replicates = replicates;
		
		for(ControlledExperiment ce : replicates){
			if(!signalSamples.contains(ce.getSignal())){
				signalSamples.add(ce.getSignal());
			}
			if(!controlSamples.contains(ce.getControl())){
				controlSamples.add(ce.getControl());
			}
			indexedReplicates.put(ce.getIndex(), ce);
			replicateIndex.put(ce, ce.getIndex());
		}
	}
	
	
	//Accessors
	
	public int getIndex(){return this.index;}
	public String getName(){return this.name;}
	public int getNumReplicates(){return this.replicates.size();}
	

}
