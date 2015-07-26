package edu.psu.compbio.seqcode.deepseq.experiments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;

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
	protected boolean estimateScaling=true;
	protected double pooledSampleControlScaling=1;
		
	public ExperimentCondition(ExptConfig c, int idx, String n, List<ControlledExperiment> reps, boolean estimateScalingFactors){
		econfig = c;
		index=idx;
		name = n;
		replicates = reps;
		numReplicates = replicates.size();
		this.estimateScaling = estimateScalingFactors;
		
		int x=0;
		for(ControlledExperiment rep : replicates){
			if(!signalSamples.contains(rep.getSignal()))
				signalSamples.add(rep.getSignal());
			if(rep.hasControl() && !controlSamples.contains(rep.getControl()))
				controlSamples.add(rep.getControl());
			replicateIndex.put(rep, rep.getIndex());
			indexedReplicate.put(rep.getIndex(), rep);
			replicateCondIndex.put(rep, x);
			
			rep.setCondition(this);
		}
		
		if(estimateScaling)
			calculateControlScalingFactors(econfig.getScalingSlidingWindow());
	}
	
	//Accessors
	public int getIndex(){return index;}
	public String getName(){return name;}
	public List<ControlledExperiment> getReplicates(){return replicates;}
	public List<Sample> getSignalSamples(){return signalSamples;}
	public List<Sample> getControlSamples(){return controlSamples;}
	public int getReplicateIndex(ControlledExperiment r){return replicateIndex.get(r);}
	public ControlledExperiment getIndexedReplicate(int id){return indexedReplicate.get(id);}
	public double getPooledSampleControlScaling(){return pooledSampleControlScaling;}
	
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
	 * Get the total weight count for signal samples from one strand
	 * @return
	 */
	public double getStrandedTotalSignalCount(char strand){
		double total=0;
		for(Sample s: signalSamples)
			total += s.getStrandedHitCount(strand);
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
	 * Get the total weight count for control samples from one strand
	 * @return
	 */
	public double getStrandedTotalControlCount(char strand){
		double total=0;
		for(Sample s: controlSamples)
			total += s.getStrandedHitCount(strand);
		return total;
	}
	
	/**
	 * Get the pooled signal fraction from all underlying signal channels
	 * @return
	 */
	public double getTotalSignalVsNoiseFrac(){
		double s=0;
		for(ControlledExperiment rep : replicates){
			s+=rep.getSignal().getHitCount()*rep.getSignalVsNoiseFraction();
		}return(s/getTotalSignalCount());
	}
	
	/**
	 * Runs the chosen scaling method for all underlying replicates (signal vs control)
	 * and for the pooled signal in this condition vs the pooled control. 
	 * 
	 * @param scalingWindowSize : size of scaling window in bp
	 */
	protected void calculateControlScalingFactors(int scalingWindowSize){
		//No point in doing any of this unless some of the replicates have controls
		if(controlSamples.size() >0){
			if(econfig.getPrintLoadingProgress())
				System.err.print("Calculating scaling factors for condition:\t"+name);
			
			ExperimentScaler scaler = new ExperimentScaler();
			Genome genome = econfig.getGenome();
			Map<Sample, List<Float>> sampleWindowCounts = new HashMap<Sample, List<Float>>();
			List<Sample> allSamples = new ArrayList<Sample>();
			allSamples.addAll(signalSamples);
			allSamples.addAll(controlSamples);
			int listSize=0;
			for(Sample samp : allSamples){
				List<Float> currSampCounts = new ArrayList<Float>();
				for(String chrom:genome.getChromList()) {
		            int chrlen = genome.getChromLength(chrom);
		            for (int start = 1; start  < chrlen - scalingWindowSize; start += scalingWindowSize) {
		                Region r = new Region(genome, chrom, start, start + scalingWindowSize);
		                currSampCounts.add(samp.countHits(r));
		            }
		        }
				sampleWindowCounts.put(samp, currSampCounts);
				listSize = currSampCounts.size();
			}
			
			//Calculate scaling factors for each replicate's signal vs control
			for(ControlledExperiment expt : getReplicates()){
				if(expt.hasControl()){
					if(econfig.getScalingBySES())
						expt.setScaling(scaler.scalingRatioBySES(sampleWindowCounts.get(expt.getSignal()), sampleWindowCounts.get(expt.getControl())));
					else if(econfig.getScalingByRegression())
						expt.setScaling(scaler.scalingRatioByRegression(sampleWindowCounts.get(expt.getSignal()), sampleWindowCounts.get(expt.getControl())));
					else if(econfig.getScalingByMedian())
						expt.setScaling(scaler.scalingRatioByMedian(sampleWindowCounts.get(expt.getSignal()), sampleWindowCounts.get(expt.getControl())));
					else
						expt.setScaling(scaler.scalingRatioByNCIS(sampleWindowCounts.get(expt.getSignal()), sampleWindowCounts.get(expt.getControl())));
				}
			}
				
			//Calculate scaling factor for pooled signal vs pooled control for this condition
			List<Float> pooledSignal = new ArrayList<Float>();
			List<Float> pooledControl = new ArrayList<Float>();
			for(int x=0; x<listSize; x++){
				float sumSig = 0;
				for(Sample s : signalSamples)
					sumSig+=sampleWindowCounts.get(s).get(x);
				pooledSignal.add(sumSig);
				
				float sumCtrl = 0;
				for(Sample s : controlSamples)
					sumCtrl+=sampleWindowCounts.get(s).get(x);
				pooledControl.add(sumCtrl);
			}
			if(econfig.getScalingBySES())
				pooledSampleControlScaling = scaler.scalingRatioBySES(pooledSignal, pooledControl);
			else if(econfig.getScalingByRegression())
				pooledSampleControlScaling = scaler.scalingRatioByRegression(pooledSignal, pooledControl);
			else
				pooledSampleControlScaling = scaler.scalingRatioByMedian(pooledSignal, pooledControl);
			
			if(econfig.getPrintLoadingProgress())
				System.err.println("\tComplete.");
		}
	}
}
