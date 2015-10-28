package edu.psu.compbio.seqcode.projects.chexmix;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;

public class CompositeModelSiteAssignment implements Comparable<CompositeModelSiteAssignment>{

	protected StrandedPoint point;
	protected Integer pointIndex;
	protected double[] totalTags; //Per-condition total tag counts for this point
	protected int[] modelComponentIndices;  //Indices of components that are active in model
	protected double[][] modelComponentResponsibilities; //Per-condition, per-component responsibility tag counts
	
	public CompositeModelSiteAssignment(StrandedPoint pt, int ptIndex, double total[], int[] modelIndices, double[][] modelResp){
		point = pt;
		pointIndex = ptIndex;
		totalTags = total;
		modelComponentIndices = modelIndices;
		modelComponentResponsibilities = modelResp;
	}
	
	public String toString(ExperimentCondition cond){
		String out = "\t"+point.getLocationString()+"\t"+totalTags[cond.getIndex()];
		for(int i=0; i<modelComponentIndices.length; i++)
			out = out+"\t"+modelComponentResponsibilities[cond.getIndex()][i];
		return out;
	}
	
	public int compareTo(CompositeModelSiteAssignment cmsa){
		return pointIndex.compareTo(cmsa.pointIndex);
	}
}
