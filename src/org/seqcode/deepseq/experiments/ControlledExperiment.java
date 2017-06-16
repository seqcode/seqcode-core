package org.seqcode.deepseq.experiments;

/**
 * ControlledExperiment combines two Samples representing signal and control
 * experiments. An object of this type may be considered as a single replicate
 * of a given ExperimentCondition. Maintains traceback links to parent
 * conditions, targets, and experiment types.
 * 
 * @author Shaun Mahony
 * @version %I%, %G%
 */
public class ControlledExperiment {

	protected ExptConfig econfig;
	protected int index;
	protected Sample signal = null; // signal channel
	protected Sample control = null; // control channel
	protected double ctrlScalingRatio = 1.0; // scaling signal channel versus
												// control channel
	protected double signalVsNoiseFraction = 0.0; // Fraction of reads assigned
													// to signal in the signal
													// channel. Init to zero to
													// parameterize Poisson
													// background model in the
													// absence of a sig/noise
													// estimate
	protected String condName;
	protected String repName;
	protected String name;
	protected ExperimentCondition myCondition = null;
	protected ExperimentTarget myTarget = null;
	protected ExperimentType myExptType = null;

	public ControlledExperiment(ExptConfig c, int idx, String cn, String rn, Sample sig, Sample ctrl) {
		econfig = c;
		index = idx;
		condName = cn;
		repName = rn;
		name = condName + ":" + repName;
		signal = sig;
		control = ctrl;

		// Initialize with read count normalization.
		// Other scaling strategies (median, regression, SES) will be
		// initialized in ExperimentCondition
		ctrlScalingRatio = control != null
				? (signal.getHitCount() / control.getHitCount()) * econfig.getFixedScalingFactor() : 1;

		// Start with assumption of no signal
		signalVsNoiseFraction = 0.0;
	}

	// Accessors
	public int getIndex() {
		return index;
	}

	public String getName() {
		return name;
	}

	public String getCondName() {
		return condName;
	}

	public String getRepName() {
		return repName;
	}

	public Sample getSignal() {
		return signal;
	}

	public Sample getControl() {
		return control;
	}

	public boolean hasControl() {
		return control != null;
	}

	public double getControlScaling() {
		return ctrlScalingRatio;
	}

	public double getSignalVsNoiseFraction() {
		return signalVsNoiseFraction;
	}

	public void setSignalVsNoiseFraction(double snf) {
		signalVsNoiseFraction = snf;
	}

	public ExperimentCondition getCondition() {
		return myCondition;
	}

	public ExperimentTarget getTarget() {
		return myTarget;
	}

	public ExperimentType getExptType() {
		return myExptType;
	}

	public void setScaling(double s) {
		ctrlScalingRatio = s;
	}

	public void setCondition(ExperimentCondition c) {
		myCondition = c;
	}

	public void setTarget(ExperimentTarget t) {
		myTarget = t;
	}

	public void setExptType(ExperimentType t) {
		myExptType = t;
	}

}
