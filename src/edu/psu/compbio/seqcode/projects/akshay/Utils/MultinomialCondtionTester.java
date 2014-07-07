package edu.psu.compbio.seqcode.projects.akshay.Utils;


import java.util.List;

import umontreal.iro.lecuyer.probdistmulti.MultinomialDist;


import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;

import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;

public class MultinomialCondtionTester {
	
	protected Config config;
	protected ExperimentSet exptSet;
	protected List<BindingEvent> features;
	protected double minFold = 2;
	protected double genomeLength;
	protected MultinomialDist mnomial;
	
	public MultinomialCondtionTester(Config con, ExperimentSet es, List<BindingEvent> features, double minFoldChange, double genomeLength ) {
		this.config = con;
		this.exptSet = es;
		this.features = features;
		this.minFold = minFoldChange;
		this.genomeLength = genomeLength;
		
		mnomial = new MultinomialDist(100, new double[2]);
		
		
	}
	
	
	public void execute(){
		
		
		for(BindingEvent cf : features){
			for(ExperimentCondition ec : exptSet.getConditions()){
				for(ControlledExperiment cr : ec.getReplicates()){
					
				}
			}
		}
	}
	
	
	
}
