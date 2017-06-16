package org.seqcode.deepseq.experiments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * ExperimentType is a collection of ControlledExperiments representing a set of
 * experiments that all share the same experimental type.
 * 
 * @author Shaun Mahony
 * @version %I%, %G%
 */
public class ExperimentType {
	protected String name;
	protected ExptConfig econfig;
	protected int index;
	protected List<ExperimentCondition> conditions = new ArrayList<ExperimentCondition>();
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected List<Sample> signalSamples = new ArrayList<Sample>();
	protected List<Sample> controlSamples = new ArrayList<Sample>();
	protected HashMap<ControlledExperiment, Integer> replicateIndex = new HashMap<ControlledExperiment, Integer>();
	protected HashMap<ControlledExperiment, Integer> replicateCondIndex = new HashMap<ControlledExperiment, Integer>();
	protected HashMap<Integer, ControlledExperiment> indexedReplicate = new HashMap<Integer, ControlledExperiment>();
	protected int numReplicates = 0;

	public ExperimentType(ExptConfig c, int idx, String n, List<ControlledExperiment> reps) {
		econfig = c;
		index = idx;
		name = n;
		replicates = reps;
		numReplicates = replicates.size();

		int x = 0;
		for (ControlledExperiment rep : replicates) {
			if (!signalSamples.contains(rep.getSignal()))
				signalSamples.add(rep.getSignal());
			if (!controlSamples.contains(rep.getControl()))
				controlSamples.add(rep.getControl());
			replicateIndex.put(rep, rep.getIndex());
			indexedReplicate.put(rep.getIndex(), rep);
			replicateCondIndex.put(rep, x);

			rep.setExptType(this);
			if (rep.getCondition() == null) {
				System.err.println("Null condition for replicate: " + rep.getName()
						+ "\nExperimentManager should initialize conditions before targets!");
				System.exit(1);
			} else {
				if (!conditions.contains(rep.getCondition()))
					conditions.add(rep.getCondition());
			}
		}
	}

	// Accessors
	public int getIndex() {
		return index;
	}

	public String getName() {
		return name;
	}

	public List<ExperimentCondition> getExptTypeConditions() {
		return conditions;
	}

	public List<ControlledExperiment> getExptTypeExperiments() {
		return replicates;
	}

	public int getReplicateIndex(ControlledExperiment r) {
		return replicateIndex.get(r);
	}

	public ControlledExperiment getIndexedReplicate(int id) {
		return indexedReplicate.get(id);
	}

}
