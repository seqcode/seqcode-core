package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import cern.jet.random.Binomial;
import cern.jet.random.ChiSquare;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.stats.Normalization;

/**
 * Test the significance of count enrichment vs control
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class EnrichmentSignificance {

	protected Config config;
	protected ExperimentSet exptSet;
	protected List<BindingEvent> features;
	protected double minFoldChange;
	protected double genomeLength;
	protected Binomial binomial;
	protected Poisson poisson;
	protected ChiSquare chisquare;
	
	public EnrichmentSignificance(Config con, ExperimentSet exptSet, List<BindingEvent> features, double minFoldChange, double genomeLength){
		this.config = con;
		this.exptSet = exptSet;
		this.features = features;
		this.minFoldChange = minFoldChange;
		this.genomeLength = genomeLength;
		binomial = new Binomial(100, .5, new DRand());
		poisson = new Poisson(1, new DRand());
		chisquare = new ChiSquare(1, new DRand());
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
		double[] repWeights = new double[exptSet.getReplicates().size()];
		for(ExperimentCondition c : exptSet.getConditions()){
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
			for(ExperimentCondition c1 : exptSet.getConditions()){
				double c1Sig = cf.getCondSigHitsFromReps(c1);
				double ctrlCountScaled = cf.getCondCtrlHitsScaledFromReps(c1);
				
				//Weighted fold difference, signal vs control
				double sigCtrlFold = 0;
				for(ControlledExperiment r : c1.getReplicates()){
					double repFold = cf.getRepCtrlHits(r)>1 ? cf.getRepSigHits(r)/(cf.getRepCtrlHits(r)*r.getControlScaling()) : cf.getRepSigHits(r);
					sigCtrlFold += repFold * repWeights[r.getIndex()];
				}
				
				//P-value, signal vs control
				double sigCtrlP = evaluateSignificance(c1Sig, ctrlCountScaled, cf.getCondTotalSigHitsFromReps(c1), c1.getMaxModelRange());
				cf.setCondSigVCtrlFold(c1, sigCtrlFold);
				cf.setCondSigVCtrlP(c1, sigCtrlP);
				
				//Log-likelihood p-value
				if(config.CALC_COMP_LL){
					cf.setLLp(c1, chisquare.cdf(cf.getLLd(c1)));
					System.out.println(String.format("%s\t%s\t%.0f\t%e",cf.getPoint().getLocationString(),c1.getName(),cf.getLLd(c1),cf.getLLp(c1)));
				}
			}
		}
		// calculate q-values, correction for multiple testing
		benjaminiHochbergCorrection(features);
	}//end of evaluateConfidence method

	
	/**
	 * Evaluate the significance using Binomial and Poisson distributions
	 */
	private double evaluateSignificance(double countA, double countB, double total, int modelWidth) {
        double pValuePoisson, pValueBalance;
		
		if(countA+countB<=0){
			return(1);
		}else{
	        try{
	
	            binomial.setNandP((int)Math.ceil(countA + countB), 1.0 / (minFoldChange + 1));
	            pValueBalance = binomial.cdf((int)Math.ceil(countB));
	
	            poisson.setMean(minFoldChange * Math.max(countB, total * (double)modelWidth / (double)genomeLength ));
	            int cA = (int)Math.ceil(countA);
	            pValuePoisson = 1 - poisson.cdf(cA) + poisson.pdf(cA);
	            
	        } catch(Exception err){
	            err.printStackTrace();
	            throw new RuntimeException(err.toString(), err);
	        }
	        return(Math.max(pValueBalance, pValuePoisson));
		}
	}
	/**
	 * Multiple hypothesis testing correction
	 */
	private void benjaminiHochbergCorrection(List<BindingEvent> features){
		double total = features.size();
		
		//Signal-vs-Control corrections by condition
		for(ExperimentCondition c : exptSet.getConditions()){
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
		
		//LL p-value corrections by condition
		if(config.CALC_COMP_LL){
			for(ExperimentCondition c : exptSet.getConditions()){
				BindingEvent.setSortingCond(c);
				Collections.sort(features, new Comparator<BindingEvent>(){
		            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareByLLPvalue(o2);}
		        });
				
				
				double rank =1.0;
				for(BindingEvent cf : features){
					cf.setCondSigVCtrlP(c, Math.min(1.0, cf.getCondSigVCtrlP(c)*(total/rank)));
					rank++;
				}
			}
		}
		
		//Finally, sort on the first condition
		BindingEvent.setSortingCond(exptSet.getConditions().get(0));
		Collections.sort(features, new Comparator<BindingEvent>(){
            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
        });
	}//end of benjaminiHochbergCorrection method
}
