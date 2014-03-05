package edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.ExperimentScaler;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.Sample;

/**
 * 
 * This is EXACTLY same as the Sample class in multigps (edu.psu.compbio.seqcode.projects.multigps.experiments), EXCEPT for
 * different Config class (Config class here is from the bayesments project and not from the multigps one)
 */
public class ControlledExperiment {

	protected Config config;
	protected int index;
	protected Sample signal=null;
	protected Sample control=null;
	protected double ctrlScalingRatio = 1.0; //scaling signal versus control
	protected double sigCount=0;   //Signal reads in the signal channel
	protected double noiseCount=0; //Noise reads in the signal channel
	protected double signalProportion= 1.0; //Fraction of reads assigned to signal
	protected BindingModel model;
	protected String condName;
	protected String repName;
	protected String feaName;
	protected String name;
	
	public ControlledExperiment(Config c, int idx, String fn, String cn, String rn, Sample sig, Sample ctrl, BindingModel initModel, boolean estimateScaling, boolean scaleByMedian){
		config = c;
		index=idx;
		condName = cn;
		repName = rn;
		feaName = fn;
		name = feaName+":"+condName+":"+repName;
		signal=sig;
		control = ctrl;
		model = initModel;
		
		if(estimateScaling){
			System.err.println("Estimating scaling ratio for "+name);
			ExperimentScaler scaler = new ExperimentScaler(signal, control);
			if(scaleByMedian)
				ctrlScalingRatio = scaler.scalingRatioByMedian(config.getScalingSlidingWindow());
			else
				ctrlScalingRatio = scaler.scalingRatioByRegression(config.getScalingSlidingWindow());
		}
	}
	
	//flushers
	
	public void flushReads(){
		this.signal.flushCounts();
		if(!(control == null))
			this.control.flushCounts();
	}
	
	//Accessors
	public int getIndex(){return index;}
	public String getName(){return name;}
	public String getFeaName(){return feaName;}
	public String getCondName(){return condName;}
	public String getRepName(){return repName;}
	public Sample getSignal(){return signal;}
	public Sample getControl(){return control;}
	public boolean hasControl(){return control!=null;}
	public double getControlScaling(){return ctrlScalingRatio;}
	public double getSigProp(){return signalProportion;}
	public BindingModel getBindingModel(){return model;}
	public void setSigNoiseCounts(double s, double n){sigCount=s; noiseCount=n; signalProportion=sigCount/(sigCount+noiseCount);}
	public double getSigCount(){return sigCount;}
	public double getNoiseCount(){return noiseCount;}
	
	// return the hits in a given regions after multiplying with the scaling ratio
	public float getSignaHitCount(Region r){
		Sample sig = this.signal;
		float count = sig.countHits(r);
		count = count*(float)this.ctrlScalingRatio;
		return count;
	}
	
	public float getSignalHitCount(int i){
		float count = this.signal.getRegionCount(i);
		count = count *(float)this.ctrlScalingRatio;
		return count;
	}
	
	public void setScaling(double s){ctrlScalingRatio = s;}
	public void setSigProp(double s){signalProportion = s;}
	public void setBindingModel(BindingModel bm){model = bm;}
		
}
