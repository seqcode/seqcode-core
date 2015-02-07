package edu.psu.compbio.seqcode.deepseq.experiments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * ExperimentCondition is a collection of ControlledExperiments representing a set of experimental replicates. 
 * Scaling replicate read counts against each other is NOT necessary under the current design philosophy.   
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentCondition {

	protected ExptConfig econfig;
	protected int index;
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected List<Sample> signalSamples = new ArrayList<Sample>();
	protected List<Sample> controlSamples = new ArrayList<Sample>();
	protected HashMap<ControlledExperiment, Integer> replicateIndex = new HashMap<ControlledExperiment, Integer>();
	protected HashMap<ControlledExperiment, Integer> replicateCondIndex = new HashMap<ControlledExperiment, Integer>();
	protected HashMap<Integer, ControlledExperiment> indexedReplicate = new HashMap<Integer, ControlledExperiment>();
	protected int numReplicates = 0;
	protected String name;
		
	public ExperimentCondition(ExptConfig c, int idx, String n, List<ControlledExperiment> reps){
		econfig = c;
		index=idx;
		name = n;
		replicates = reps;
		numReplicates = replicates.size();
		
		int x=0;
		for(ControlledExperiment rep : replicates){
			if(!signalSamples.contains(rep.getSignal()))
				signalSamples.add(rep.getSignal());
			if(!controlSamples.contains(rep.getControl()))
				controlSamples.add(rep.getControl());
			replicateIndex.put(rep, rep.getIndex());
			indexedReplicate.put(rep.getIndex(), rep);
			replicateCondIndex.put(rep, x);
			
			rep.setCondition(this);
		}
	}
	
	//Accessors
	public int getIndex(){return index;}
	public String getName(){return name;}
	public List<ControlledExperiment> getReplicates(){return replicates;}
	public List<Sample> getSignalSamples(){return signalSamples;}
	public List<Sample> getControlSamples(){return controlSamples;}
	public int getReplicateIndex(ControlledExperiment r){return replicateIndex.get(r);}
	public ControlledExperiment getIndexedReplicate(int id){return indexedReplicate.get(id);}
	
	
	/**
	 * Get the total weight count for signal samples
	 * @return
	 */
	public double getTotalSignalCount(){
		double total=0;
		for(Sample s: signalSamples)
			total += s.getHitCount();
		return total;
	}
	
	/**
	 * Get the total weight count for the estimated background components in signal samples
	 * @return
	 */
	public double getTotalSignalEstBackCount(){
		double total=0;
		for(ControlledExperiment r : replicates){
			total+=r.getNoiseCount();
		}
		return total;
	}
	
	/**
	 * Get the total weight count for signal samples
	 * @return
	 */
	public double getTotalControlCount(){
		double total=0;
		for(Sample s: controlSamples)
			total += s.getHitCount();
		return total;
	}
	/**
	 * Get the fraction of reads from all replicates in signal
	 * @return
	 */
	public double getSigProp(){
		double s=0, n=0;
		for(ControlledExperiment rep : replicates){
			s+=rep.getSigCount();
			n+=rep.getNoiseCount();
		}return(s/(s+n));
	}
	
}
