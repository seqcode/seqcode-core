package edu.psu.compbio.seqcode.projects.multigps.experiments;

import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.ExperimentScaler;

/**
 * ControlledExperiment combines two Samples representing signal and control experiments. 
 * An object of this type may be considered as a single replicate of a given Condition.
 * This class contains a replicate-specific BindingModel.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
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
	protected String name;
	
	public ControlledExperiment(Config c, int idx, String cn, String rn, Sample sig, Sample ctrl, BindingModel initModel, boolean estimateScaling, boolean scaleByMedian){
		config = c;
		index=idx;
		condName = cn;
		repName = rn;
		name = condName+":"+repName;
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
	
	//Accessors
	public int getIndex(){return index;}
	public String getName(){return name;}
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
	
	public void setScaling(double s){ctrlScalingRatio = s;}
	public void setSigProp(double s){signalProportion = s;}
	public void setBindingModel(BindingModel bm){model = bm;}
		
}
