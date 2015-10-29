package edu.psu.compbio.seqcode.projects.chexmix;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;

public class CompositeModelSiteAssignment{

	protected StrandedPoint point;
	protected double[] totalTags; //Per-condition total tag counts for this point
	protected int[] modelComponentIndices;  //Indices of components that are active in model
	protected double[][] modelComponentResponsibilities; //Per-condition, per-component responsibility tag counts
	
	public CompositeModelSiteAssignment(StrandedPoint pt, double total[], int[] modelIndices, double[][] modelResp){
		point = pt;
		totalTags = total;
		modelComponentIndices = modelIndices;
		modelComponentResponsibilities = modelResp;
	}
	
	public String toString(ExperimentCondition cond){
		String out = point.getLocationString()+String.format("\t%.2f", totalTags[cond.getIndex()]);
		for(int i=0; i<modelComponentIndices.length; i++)
			out = out+String.format("\t%.2f", modelComponentResponsibilities[cond.getIndex()][i]);
		return out;
	}
	
}
