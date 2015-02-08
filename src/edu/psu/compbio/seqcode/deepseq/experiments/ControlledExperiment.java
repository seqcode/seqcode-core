package edu.psu.compbio.seqcode.deepseq.experiments;


/**
 * ControlledExperiment combines two Samples representing signal and control experiments. 
 * An object of this type may be considered as a single replicate of a given ExperimentCondition.
 * Maintains traceback links to parent conditions, targets, and experiment types.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ControlledExperiment {

	protected ExptConfig econfig;
	protected int index;
	protected Sample signal=null;
	protected Sample control=null;
	protected double ctrlScalingRatio = 1.0; //scaling signal versus control
	protected double sigCount=0;   //Signal reads in the signal channel
	protected double noiseCount=0; //Noise reads in the signal channel
	protected double signalProportion= 0.0; //Fraction of reads assigned to signal. Init to zero to parameterize Poisson background model in the absence of a sig/noise estimate
	protected boolean estimateScaling=true;
	protected String condName;
	protected String repName;
	protected String name;
	protected ExperimentCondition myCondition=null;
	protected ExperimentTarget myTarget = null;
	protected ExperimentType myExptType= null;
	
	public ControlledExperiment(ExptConfig c, int idx, String cn, String rn, Sample sig, Sample ctrl, boolean estimateScaling){
		econfig = c;
		index=idx;
		condName = cn;
		repName = rn;
		name = condName+":"+repName;
		signal=sig;
		control = ctrl;
		this.estimateScaling = estimateScaling;
		
		if(estimateScaling){
			ExperimentScaler scaler = new ExperimentScaler(signal, control);
			if(econfig.getScalingBySES())
				ctrlScalingRatio = scaler.scalingRatioBySES(econfig.getScalingSlidingWindow());
			else if(econfig.getScalingByRegression())
				ctrlScalingRatio = scaler.scalingRatioByRegression(econfig.getScalingSlidingWindow());
			else
				ctrlScalingRatio = scaler.scalingRatioByMedian(econfig.getScalingSlidingWindow());
			
			//I don't trust that this works all the time
			//signalProportion = 1-scaler.calculateBackgroundFromScalingRatio();
			//sigCount = signalProportion*signal.getHitCount();
			//noiseCount = (1-signalProportion)*signal.getHitCount();
		}
		//Start with assumption of no signal
		sigCount = 0.0;
		noiseCount = signal.getHitCount();
		signalProportion = 0.0;
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
	public void setSigNoiseCounts(double s, double n){sigCount=s; noiseCount=n; signalProportion=sigCount/(sigCount+noiseCount);}
	public double getSigCount(){return sigCount;}
	public double getNoiseCount(){return noiseCount;}
	public ExperimentCondition getCondition(){return myCondition;}
	public ExperimentTarget getTarget(){return myTarget;}
	public ExperimentType getExptType(){return myExptType;}
	
	public void setScaling(double s){ctrlScalingRatio = s;}
	public void setCondition(ExperimentCondition c){myCondition = c;}
	public void setTarget(ExperimentTarget t){myTarget = t;}
	public void setExptType(ExperimentType t){myExptType = t;}
	
	/**
	 * Call linear count correction method in signal and recalculate necessary variables here 
	 * @param perBaseScaling
	 */
	public void correctSignalCounts(float perBaseScaling){
		signal.linearCountCorrection(perBaseScaling);
		if(estimateScaling){
			ExperimentScaler scaler = new ExperimentScaler(signal, control);
			if(econfig.getScalingBySES())
				ctrlScalingRatio = scaler.scalingRatioBySES(econfig.getScalingSlidingWindow());
			else if(econfig.getScalingByRegression())
				ctrlScalingRatio = scaler.scalingRatioByRegression(econfig.getScalingSlidingWindow());
			else
				ctrlScalingRatio = scaler.scalingRatioByMedian(econfig.getScalingSlidingWindow());
			
			signalProportion = 1-scaler.calculateBackgroundFromScalingRatio();
			sigCount = signalProportion*signal.getHitCount();
			noiseCount = (1-signalProportion)*signal.getHitCount();
		}
	}
		
}
