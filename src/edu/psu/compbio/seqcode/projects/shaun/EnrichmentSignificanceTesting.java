package edu.psu.compbio.seqcode.projects.shaun;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;

/**
 * Test the significance of count enrichment vs control
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class EnrichmentSignificanceTesting {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager exptMan;
	protected List<BindingEvent> features;
	protected double minFoldChange;
	protected double genomeLength;
	protected Binomial binomial;
	
	public EnrichmentSignificanceTesting(GenomeConfig gcon, ExptConfig econ, List<BindingEvent> features, double minFoldChange, double genomeLength){
		this.gconfig = gcon;
		this.econfig = econ;
		this.exptMan = new ExperimentManager(econfig);
		this.features = features;
		this.minFoldChange = minFoldChange;
		this.genomeLength = genomeLength;
		binomial = new Binomial(100, .5, new DRand());
	}

	/**
	 * Evaluate the significance of a set of EnrichedFeatures using Binomial and Poisson tests
	 * Assumes that the counts in the features are not already scaled
	 * @param features List of EnrichedFeatures
	 * @param minFoldChange
	 * @param modelWidth
	 * @param totalSignalCount
	 * @param genomeLength
	 */
	public void execute() {

		//Calculate relative replicate weights
		double[] repWeights = new double[exptMan.getReplicates().size()];
		for(ExperimentCondition c : exptMan.getConditions()){
			double totalSig =0;
			for(ControlledExperiment r : c.getReplicates())
				totalSig += r.getSigCount();
			for(ControlledExperiment r : c.getReplicates())
				repWeights[r.getIndex()]=r.getSigCount()/totalSig;
		}
		
		//Compute fold difference, log-likelihood, and p-values. 
		//To be clear about what each variable is measuring:
		// sigCtrlFold: weighted average of per-replicate fold difference between signal reads assigned to event and scaled control reads assigned to event.
		// sigCtrlP: Outcome of binomial test between sum of per-replicate signal read counts at event and sum of scaled per-replicate control counts at event. If multiple replicates use the same control, the counts are used redundantly since this is equivalent to summing the scaling factors across replicates. 
		// LL: log-likelihood loss if component was eliminated from model.
		// LLp: Chi-square distributed p-value corresponding to LL.  
		for (BindingEvent cf: features){
			for(ExperimentCondition c1 : exptMan.getConditions()){
				double c1Sig = cf.getCondSigHitsFromReps(c1);
				double ctrlCountScaled = cf.getCondCtrlHitsScaledFromReps(c1);
				
				//Weighted fold difference, signal vs control
				double sigCtrlFold = c1Sig / ctrlCountScaled;
				
				//P-value, signal vs control
				double sigCtrlP = evaluateSignificance(c1Sig, ctrlCountScaled);
				cf.setCondSigVCtrlFold(c1, sigCtrlFold);
				cf.setCondSigVCtrlP(c1, sigCtrlP);
				//System.out.println(c1.getName()+"\t"+c1Sig+"\t"+ctrlCountScaled+"\t"+sigCtrlFold+"\t"+sigCtrlP);
			}
		}
		// calculate q-values, correction for multiple testing
		benjaminiHochbergCorrection(features);
	}//end of evaluateConfidence method

	
	/**
	 * Evaluate the significance using Binomial and Poisson distributions
	 */
	private double evaluateSignificance(double countA, double countB) {
        double pValuePoisson, pValueBalance;
		
		if(countA+countB<=0){
			return(1);
		}else{
	        try{
	
	            binomial.setNandP((int)Math.ceil(countA + countB), 1.0 / (minFoldChange + 1));
	            pValueBalance = binomial.cdf((int)Math.ceil(countB));
	
	            
	        } catch(Exception err){
	            err.printStackTrace();
	            throw new RuntimeException(err.toString(), err);
	        }
	        return(pValueBalance);
		}
	}
	/**
	 * Multiple hypothesis testing correction
	 */
	private void benjaminiHochbergCorrection(List<BindingEvent> features){
		double total = features.size();
		
		//Signal-vs-Control corrections by condition
		for(ExperimentCondition c : exptMan.getConditions()){
			BindingEvent.setSortingCond(c);
			Collections.sort(features, new Comparator<BindingEvent>(){
	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
	        });
			
			
			double rank =1.0;
			for(BindingEvent cf : features){
				cf.setCondSigVCtrlP(c, Math.min(1.0, cf.getCondSigVCtrlP(c)*(total/rank)));
				rank++;
			}
		}
		
		//Finally, sort on the first condition
		BindingEvent.setSortingCond(exptMan.getConditions().get(0));
		Collections.sort(features, new Comparator<BindingEvent>(){
            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
        });
	}//end of benjaminiHochbergCorrection method
}
