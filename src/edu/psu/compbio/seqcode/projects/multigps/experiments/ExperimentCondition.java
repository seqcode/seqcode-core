package edu.psu.compbio.seqcode.projects.multigps.experiments;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

/**
 * ExperimentCondition is a collection of ControlledExperiments representing a set of experimental replicates. 
 * Scaling replicate read counts against each other is NOT necessary under the current design philosophy.   
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentCondition {

	protected Config config;
	protected int index;
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected List<Sample> signalSamples = new ArrayList<Sample>();
	protected List<Sample> controlSamples = new ArrayList<Sample>();
	protected List<BindingEvent> events;
	protected HashMap<ControlledExperiment, Integer> replicateIndex = new HashMap<ControlledExperiment, Integer>();
	protected HashMap<ControlledExperiment, Integer> replicateCondIndex = new HashMap<ControlledExperiment, Integer>();
	protected HashMap<Integer, ControlledExperiment> indexedReplicate = new HashMap<Integer, ControlledExperiment>();
	protected WeightMatrix motif = null; //Motif for motif priors
	protected WeightMatrix freqMatrix = null; //Freq matrix version of motif for motif priors (easier alignment)
	protected int motifOffset=0; //Offset in motif alignment (against other conditions)
	protected int numReplicates = 0;
	protected int maxInfluenceRange=0; //Maximum influence range of replicate binding models 
	protected String name;
		
	public ExperimentCondition(Config c, int idx, String n, List<ControlledExperiment> reps){
		config = c;
		index=idx;
		name = n;
		replicates = reps;
		numReplicates = replicates.size();
		events = new ArrayList<BindingEvent>();
		
		int x=0;
		for(ControlledExperiment rep : replicates){
			if(!signalSamples.contains(rep.getSignal()))
				signalSamples.add(rep.getSignal());
			if(!controlSamples.contains(rep.getControl()))
				controlSamples.add(rep.getControl());
			replicateIndex.put(rep, rep.getIndex());
			indexedReplicate.put(rep.getIndex(), rep);
			replicateCondIndex.put(rep, x);
		}
		this.updateMaxInfluenceRange();
	}
	
	//Accessors
	public int getIndex(){return index;}
	public String getName(){return name;}
	public List<ControlledExperiment> getReplicates(){return replicates;}
	public int getReplicateIndex(ControlledExperiment r){return replicateIndex.get(r);}
	public ControlledExperiment getIndexedReplicate(int id){return indexedReplicate.get(id);}
	public List<BindingEvent> getBindingEvents(){return events;}
	public WeightMatrix getMotif(){return motif;}
	public WeightMatrix getFreqMatrix(){return freqMatrix;}
	public int getMotifOffset(){return motifOffset;}
	public int getMaxInfluenceRange(){return maxInfluenceRange;}
	
	//Setters
	public void setBindingEvents(List<BindingEvent> e){events=e;}
	public void addBindingEvent(BindingEvent e){events.add(e);}
	public void setMotif(WeightMatrix m){motif = m;}
	public void setFreqMatrix(WeightMatrix m){freqMatrix = m;}
	public void setMotifOffset(int o){motifOffset=o;}
	public void updateMaxInfluenceRange(){
		maxInfluenceRange=0; 
		for(ControlledExperiment rep : replicates){
			if(rep.getBindingModel().getInfluenceRange()>maxInfluenceRange)
				maxInfluenceRange=rep.getBindingModel().getInfluenceRange();
		}
	}
	
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
	/**
	 * Get the maximum model width
	 * @return
	 */
	public int getMaxModelRange(){
		int max=0;
		for(ControlledExperiment ce : replicates)
			if(ce.getBindingModel().getInfluenceRange()>max)
				max = ce.getBindingModel().getInfluenceRange();
		return max;
	}
	
	/**
	 * Print the events
	 * @param outRoot
	 */
	public void printEvents(String outRoot){
		try{
			String outName = outRoot+"."+name+".events";
			FileWriter fw = new FileWriter(outName);
			
			fw.write(BindingEvent.fullHeadString()+"\n");
			for(BindingEvent e : events){
				fw.write(e.toString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Print the replicate counts at events
	 * @param outRoot
	 */
	public void printReplicateCounts(String outRoot){
		try{
			String outName = outRoot+"."+name+".repcounts";
			FileWriter fw = new FileWriter(outName);
			
			fw.write(BindingEvent.repCountHeadString()+"\n");
			for(BindingEvent e : events){
				fw.write(e.getRepCountString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
