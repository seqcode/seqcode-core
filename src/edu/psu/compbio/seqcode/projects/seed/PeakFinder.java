package edu.psu.compbio.seqcode.projects.seed;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.projects.seed.SEEDConfig.PeakFindingMethod;
import edu.psu.compbio.seqcode.projects.seed.features.EnrichedFeature;
import edu.psu.compbio.seqcode.projects.seed.features.EnrichedPeakFeature;
import edu.psu.compbio.seqcode.projects.seed.features.Feature;

public class PeakFinder extends DomainFinder {
	public static String version = "0.1";

	public PeakFinder(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, ExperimentManager man) {
		super(gcon, econ, scon, man);
	}
	
	/**
	 * Return the class name
	 */
	public String getProgramName(){
		return "edu.psu.compbio.seqcode.projects.seed.PeakFinder";
	}
	
	/**
	 * Return a thread for this implementation
	 */
	public FeatureDetectionThread getMyThread(List<Region> regs){
		return new PeakFinderThread(regs);
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		System.err.println("PeakFinder version "+PeakFinder.version+"\n\n");
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		SEEDConfig scon = new SEEDConfig(gcon, args);
		
		if(scon.helpWanted()){
			//System.out.println(PeakFinder.getDomainFinderArgs());
			System.err.println(gcon.getArgsList()+
					econ.getArgsList()+
					scon.getArgsList());
		}else{
			ExperimentManager man = new ExperimentManager(econ);
			PeakFinder finder = new PeakFinder(gcon, econ, scon, man);
			System.err.println("\nBeginning peak finding...");
			finder.execute();
			man.close();
		}
	}
	
	/**
	 * Return a help string
	 * @return
	 */
	public static String getPeakFinderArgs() {
		//TODO
		return("TODO");
	}
	
	/**
	 * Multiple-hypothesis correction. Print the output files. 
	 * Assumes that all peaks have been found.
	 * @return : final features
	 */
	public Map<ExperimentCondition, List<Feature>> postProcess() {
		System.err.println("\nPeak finding complete.");
		for(ExperimentCondition cond : manager.getConditions()){
			stats.benjaminiHochbergCorrection(features.get(cond));
		}
		
		Map<ExperimentCondition, List<Feature>> signifFeatures = this.filter(features, sconfig.perBinBinomialPThres, true);
       	
		//All domains
		this.printEventsFile(features, ".all.peaks");
        
		//Filtered by q-value
		this.printEventsFile(signifFeatures, ".p"+sconfig.perBinBinomialPThres+".peaks");
		
		//Summarize
		for(ExperimentCondition cond : manager.getConditions())
			System.err.println(cond.getName()+"\t"+features.get(cond).size()+" peaks\t"+signifFeatures.get(cond).size()+" peaks below threshold.");
		
		return features;
	}
	
	
	/**
	 * PeakFinderThread: thread that searches for domains
	 * @author mahony
	 *
	 */
	public class PeakFinderThread extends DomainFinderThread {

		public PeakFinderThread(List<Region> regs) {
			super(regs);
		}
		
		/**
		 * findFeatures: find potentially enriched domains on the genome by comparing tag counts to a background model,
		 * and find the peaks of enrichment within these domains using a choice of methods. 
		 * 
		 * @param subRegion : region to run analysis on
		 * @return Map of Lists of Features in the subRegion, Indexed by ExperimentCondition 
		 */
		public Map<ExperimentCondition, List<Feature>> findFeatures(Region subRegion) {
			//Get the enriched domains using the functionality from the superclass (DomainFinderThread) 
			//and the over-ridden processDomain methods below 
			Map<ExperimentCondition, List<Feature>> peaks = super.findFeatures(subRegion);
			return peaks;
		}
		
		/**
		 * ProcessDomains: 
		 *  - Trims feature coordinates back to agree with overlapping hits. 
		 *  - Counts hits in each feature, per sample 
		 *  - Finds the peak position in the domain according to a choice of methods.
		 *  
		 * Specifically, the procedure characterizes enriched domains using the functionality in the superclass (DomainFinder),
		 * and then finds peaks using either:
		 * 		- Point of maximum density
		 * 		- Point of balance between positive and negative strand densities
		 * 		- Point of maximum likelihood, given an empirical read distribution (stranded) 
		 * 
		 * @param currFeatures
		 * @param current region
		 * @return : Lists of EnrichedFeatures, indexed by condition
		 */
		protected Map<ExperimentCondition, List<EnrichedFeature>> processDomains(Map<ExperimentCondition,List<EnrichedFeature>> currFeatures, Region currSubRegion){
			Map<ExperimentCondition, List<EnrichedFeature>> peakFeatures = new HashMap<ExperimentCondition, List<EnrichedFeature>>();
			for(ExperimentCondition cond : manager.getConditions())
				peakFeatures.put(cond, new ArrayList<EnrichedFeature>());
			
			for(ExperimentCondition currCondition : manager.getConditions()){
				for(EnrichedFeature currDomain : currFeatures.get(currCondition)){
					Map<Sample, List<StrandedBaseCount>> fHitsPos = overlappingHits(hitsPos, currDomain);
					Map<Sample, List<StrandedBaseCount>> fHitsNeg = overlappingHits(hitsNeg, currDomain);
					
					//Trim the coordinates
					trimFeature(currDomain, fHitsPos, fHitsNeg, currCondition);

					//Quantify the feature in each Sample and in the condition in which it was found
					quantifyFeature(currDomain, fHitsPos, fHitsNeg, currCondition);
					
					//Find the peaks
					EnrichedPeakFeature peak;
					if(sconfig.getPeakFindingApproach() == PeakFindingMethod.LRBALANCE)
	                    peak = findPeakLRBalance(fHitsPos, fHitsNeg, (EnrichedFeature)currDomain, currCondition);
	                //else if(sconfig.getPeakFindingApproach() == PeakFindingMethod.TAGDISTRIBUTION && sconfig.getEventModel()!=null)
	                //    peak = findPeakWithBindingModel(fHitsPos, fHitsNeg, currDomain, currCondition);
	                else if(sconfig.getPeakFindingApproach() == PeakFindingMethod.MAXDENSITY)
	                    peak = findPeakMaxHit(fHitsPos, fHitsNeg, (EnrichedFeature)currDomain, currCondition);
	                else{
	                	System.err.println("PeakFindingMethod "+sconfig.getPeakFindingApproach()+" not defined. Using max density.");
	                	peak = findPeakMaxHit(fHitsPos, fHitsNeg, (EnrichedFeature)currDomain, currCondition);
	                }
					if(peak!=null)
						peakFeatures.get(currCondition).add(peak);
				}
			}
			return(peakFeatures);
		}
		
		/**
		 * Find the peak locations based on maximum overlapping read counts. 
	     * 
		 * @return : Feature (EnrichedPeakFeature)
		 */
		protected EnrichedPeakFeature findPeakMaxHit(Map<Sample, List<StrandedBaseCount>> fHitsPos, Map<Sample, List<StrandedBaseCount>> fHitsNeg, 
									EnrichedFeature domain, ExperimentCondition currCondition){
			float [] sum = new float[domain.getCoords().getWidth()+1];
			for(int s=0; s<=domain.getCoords().getWidth(); s++){sum[s]=0; }

			for(Sample s : currCondition.getSignalSamples()){
				if(!strandedEventDetection || domain.getCoords().getStrand()=='+')
					for(StrandedBaseCount h : fHitsPos.get(s)){
						int start = getLeft(h)-domain.getCoords().getStart(); 
			            int stop= getRight(h)-domain.getCoords().getStart();
			            for(int i=start; i<=stop; i++)
			                if(i>=0 && i<sum.length)
			                    sum[i]+=h.getCount();
					}
				if(!strandedEventDetection || domain.getCoords().getStrand()=='-')
					for(StrandedBaseCount h : fHitsNeg.get(s)){
						int start = getLeft(h)-domain.getCoords().getStart(); 
			            int stop= getRight(h)-domain.getCoords().getStart();
			            for(int i=start; i<=stop; i++)
			                if(i>=0 && i<sum.length)
			                    sum[i]+=h.getCount();
					}
			}
			float max = 0; int maxPos = -1;
			for(int s=0; s<sum.length; s++){
				if(sum[s]>max){
					max= sum[s];
					maxPos=s;
				}
			}
			Point p = new Point(gen, domain.getCoords().getChrom(), maxPos+domain.getCoords().getStart());
			EnrichedPeakFeature epf=null;
			try {
				epf = new EnrichedPeakFeature(domain.getCoords(), p, domain.getSampleCountsPos(), domain.getSampleCountsNeg(), 
																	domain.getSignalCount(), domain.getControlCount(), domain.getScore(), max);
			} catch (Exception e) {}
			return(epf);
		}
		
		/**
		 * Find the peak locations based on left-right balance of forward/reverse reads. 
		 * Motivated by SISSRS.
		 * 
		 * @param hits
		 * @param coords
		 * @return
		 */
		protected EnrichedPeakFeature findPeakLRBalance(Map<Sample, List<StrandedBaseCount>> fHitsPos, Map<Sample, List<StrandedBaseCount>> fHitsNeg, 
									EnrichedFeature domain, ExperimentCondition currCondition){
			float [] forward = new float [domain.getCoords().getWidth()+1];
			float [] reverse = new float [domain.getCoords().getWidth()+1];
			float [] density = new float[domain.getCoords().getWidth()+1];
			for(int s=0; s<=domain.getCoords().getWidth(); s++){forward[s]=0; reverse[s]=0;}
			
			if(strandedEventDetection){
				System.err.println("Don't use this method for single-stranded event detection!");
				System.exit(1);
			}
			
			for(Sample s : currCondition.getSignalSamples()){
				for(StrandedBaseCount h : fHitsPos.get(s)){
					int fivePrime = getShifted5Prime(h)-domain.getCoords().getStart();
					if(fivePrime>=0 && fivePrime<forward.length)
						forward[fivePrime]+=h.getCount();
					int start = getLeft(h)-domain.getCoords().getStart(); 
		            int stop= getRight(h)-domain.getCoords().getStart();
		            for(int i=start; i<stop; i++)
		                if(i>=0 && i<density.length)
		                    density[i]+=h.getCount();
				}
				for(StrandedBaseCount h : fHitsNeg.get(s)){
					int fivePrime = getShifted5Prime(h)-domain.getCoords().getStart();
					if(fivePrime>=0 && fivePrime<reverse.length)
						reverse[fivePrime]+=h.getCount();
					int start = getLeft(h)-domain.getCoords().getStart(); 
		            int stop= getRight(h)-domain.getCoords().getStart();
		            for(int i=start; i<stop; i++)
		                if(i>=0 && i<density.length)
		                    density[i]+=h.getCount();
				}
			}
			
			int minBal = Integer.MAX_VALUE, minPos = -1;
			for(int s=0; s<=domain.getCoords().getWidth(); s++){
				int left=0, right=0;
				for(int l=0; l<=s; l++){left+=forward[l];}
				for(int r=domain.getCoords().getWidth(); r>s; r--){right+=reverse[r];}
				if(Math.abs(left-right)<minBal){
					minBal =Math.abs(left-right);
					minPos = s;
				}
			}
			Point p = new Point(gen, domain.getCoords().getChrom(), minPos+domain.getCoords().getStart());
			EnrichedPeakFeature epf=null;
			try {
				epf = new EnrichedPeakFeature(domain.getCoords(), p, domain.getSampleCountsPos(), domain.getSampleCountsNeg(), 
																	domain.getSignalCount(), domain.getControlCount(), 
																	domain.getScore(), density[minPos]);
			} catch (Exception e) {}
			return(epf);
		}
		
		
		/* Find exact peak using a BindingModel */
/*		protected EnrichedFeature findPeakWithBindingModel(Map<Sample, List<StrandedBaseCount>> fHitsPos, Map<Sample, List<StrandedBaseCount>> fHitsNeg, 
								EnrichedFeature domain, ExperimentCondition currCondition, BindingModel model){
			float [] p = new float[domain.getCoords().getWidth()+1];
			for(int s=0; s<=domain.getCoords().getWidth(); s++){p[s]=0; }
			
			int maxPos=0; float maxScore=0;
			//Populate the probability landscape - need to correct for strandedness of events here.
			for(Sample s : currCondition.getSignalSamples()){
				for(StrandedBaseCount h : fHitsPos.get(s)){
					int fivePrime = getShifted5Prime(h)-domain.getCoords().getStart();
					if(fivePrime>=0 && fivePrime<p.length){
						for(int i=Math.max(model.getMin()+fivePrime, 0); i<=Math.min(domain.getCoords().getWidth(), fivePrime+model.getMax()); i++)
							p[i]+=model.probability(i-fivePrime);
					}
				}
				for(StrandedBaseCount h : fHitsNeg.get(s)){
					int fivePrime = getShifted5Prime(h)-domain.getCoords().getStart();
					if(fivePrime>=0 && fivePrime<p.length){
						for(int i=Math.max(fivePrime-model.getMax(), 0); i<=Math.min(domain.getCoords().getWidth(), fivePrime-model.getMin()); i++)
							p[i]+=model.probability(fivePrime-i);
					}
				}
			}
			for(int k=0; k<=domain.getCoords().getWidth(); k++)
				if(p[k]>maxScore){maxScore=p[k]; maxPos=k;}
				
			Point pt = new Point(gen, domain.getCoords().getChrom(), domain.getCoords().getStart()+maxPos);
			EnrichedPeakFeature epf=null;
			try {
				epf = new EnrichedPeakFeature(domain.getCoords(), pt, domain.getSampleCountsPos(), domain.getSampleCountsNeg(), 
																	domain.getSignalCount(), domain.getControlCount(), 
																	domain.getScore(), maxScore);
			} catch (Exception e) {}
			return(epf);
		}
		*/
	}

	
}
