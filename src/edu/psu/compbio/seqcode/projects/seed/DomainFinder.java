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
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.projects.seed.features.EnrichedFeature;
import edu.psu.compbio.seqcode.projects.seed.features.Feature;
import edu.psu.compbio.seqcode.projects.seed.stats.FeatureStatistics;

/**
 * DomainFinder: A FeatureDetection implementation that scans the genome for statistically enriched blocks (domains)
 * 
 * @author mahony
 *
 */
public class DomainFinder extends FeatureDetection {

	public static String version = "0.1";
	protected FeatureStatistics stats;
	
	public DomainFinder(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, ExperimentManager man) {
		super(gcon, econ, scon, man);
		stats = new FeatureStatistics();
	}

	/**
	 * Return the class name
	 */
	public String getProgramName(){
		return "edu.psu.compbio.seqcode.projects.seed.DomainFinder";
	}
	
	/**
	 * Return a thread for this implementation
	 */
	public FeatureDetectionThread getMyThread(List<Region> regs){
		return new DomainFinderThread(regs);
	}
	
	/**
	 * Main command-line executable method
	 * @param args
	 */
	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		System.err.println("DomainFinder version "+DomainFinder.version+"\n\n");
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		SEEDConfig scon = new SEEDConfig(gcon, args);
		
		if(scon.helpWanted()){
			//System.out.println(DomainFinder.getDomainFinderArgs());
			System.err.println(gcon.getArgsList()+
					econ.getArgsList()+
					scon.getArgsList());
		}else{
			ExperimentManager man = new ExperimentManager(econ);
			DomainFinder finder = new DomainFinder(gcon, econ, scon, man);
			System.err.println("\nBeginning domain finding...");
			finder.execute();
			man.close();
		}
	}
	
	/**
	 * Return a help string
	 * @return
	 */
	public static String getDomainFinderArgs() {
		//TODO
		return("TODO");
	}
	
	/**
	 * Multiple-hypothesis correction. Print the output files. 
	 * Assumes that all domains have been found.
	 * @return : final features
	 */
	public Map<ExperimentCondition, List<Feature>> postProcess() {
		System.err.println("\nDomain finding complete.");
		for(ExperimentCondition cond : manager.getConditions()){
			stats.benjaminiHochbergCorrection(features.get(cond));
		}
		
		Map<ExperimentCondition, List<Feature>> signifFeatures = this.filter(features, sconfig.perBinBinomialPThres, true);
       	
		//All domains
		this.printEventsFile(features, ".all.domains");
        
		//Filtered by q-value
		this.printEventsFile(signifFeatures, ".p"+sconfig.perBinBinomialPThres+".domains");
		
		//Summarize
		for(ExperimentCondition cond : manager.getConditions())
			System.err.println(cond.getName()+"\t"+features.get(cond).size()+" domains\t"+signifFeatures.get(cond).size()+" domains below threshold.");
		
		return features;
	}

	
	/**
	 * DomainFinderThread: thread that searches for domains
	 * @author mahony
	 *
	 */
	public class DomainFinderThread extends FeatureDetectionThread {

		public DomainFinderThread(List<Region> regs) {
			super(regs);
		}

		/**
		 * findFeatures: find potentially enriched domains on the genome by comparing tag counts to a background model
		 * 
		 * Specifically, the procedure requires that the signal count in a bin passes all PoissonBackground thresholds, 
		 * and the signal count vs control count (if exists) Binomial test passes an uncorrected p-value threshold.  
		 * Neighboring bins that pass the thresholds are merged.
		 * 
		 * This method can assume that the hits and landscape data structures have been initialized for a given region 
		 * - for each Sample these data structures contain, respectively: the StrandedBaseCounts hits (sorted), and the 
		 * binned tag density (after any shifting and extending).
		 * 
		 * @param subRegion : region to run analysis on
		 * @return Map of Lists of Features in the subRegion, Indexed by ExperimentCondition 
		 */
		public Map<ExperimentCondition, List<Feature>> findFeatures(Region subRegion) {
			Map<ExperimentCondition, List<Feature>> results = new HashMap<ExperimentCondition, List<Feature>>();
			for(ExperimentCondition cond : manager.getConditions())
				results.put(cond, new ArrayList<Feature>());
			Map<ExperimentCondition, List<EnrichedFeature>> currFeatures = new HashMap<ExperimentCondition, List<EnrichedFeature>>();
			for(ExperimentCondition cond : manager.getConditions())
				currFeatures.put(cond, new ArrayList<EnrichedFeature>());
			
			for(ExperimentCondition cond : manager.getConditions()){
				int numStrandIter = strandedEventDetection ? 2 : 1;
				for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
					EnrichedFeature lastFeature=null;
	                //If stranded peak-finding, run over both strands separately
	                char str = !strandedEventDetection ? '.' : (stranditer==1 ? '+' : '-');
	                
	                float [] condSigCounts = getConditionCounts(cond, landscape, str, true);
	                float [] condCtrlCounts = cond.getControlSamples().size()>0 ?
	                							null : getConditionCounts(cond, landscape, str, false);
	                
	                //Scan regions
	                int currBin=0;
	                for(int i=subRegion.getStart(); i<subRegion.getEnd()-sconfig.getBinWidth(); i+=sconfig.getBinStep()){
	                	float sigCounts=condSigCounts[currBin];
	                	//If there's no control, use the expected bin count given random distribution
	                	float ctrlCounts= condCtrlCounts==null ? conditionBackgrounds.get(cond).getGenomicExpectedCount()
	                			: condCtrlCounts[currBin]*(float)(cond.getPooledSampleControlScaling());
	                	
                        //First Test: is the read count above the genome-wide thresholds? 
                        if(conditionBackgrounds.get(cond).passesGenomicThreshold((int)sigCounts, str)){
                        	//Second Test: refresh all thresholds & test again
                    		conditionBackgrounds.get(cond).updateModels(subRegion, i-subRegion.getStart(), condSigCounts, condCtrlCounts, sconfig.getBinStep());
                        	if(conditionBackgrounds.get(cond).passesAllThresholds((int)sigCounts, str)){
                        		//Third Test: Binomial test between signal & (scaled) control 
                        		double pval = stats.binomialPValue(ctrlCounts, (sigCounts+ctrlCounts), sconfig.getMinSigCtrlFoldDifference());
                        		if(pval < sconfig.getPerBinBinomialPThres()){
                        			//Add new event or append to the last one.
                        			Region currWin = strandedEventDetection ? 
                        					new StrandedRegion(gen, subRegion.getChrom(), i, i+sconfig.getBinWidth()-1, str) :
                        					new Region(gen, subRegion.getChrom(), i, i+sconfig.getBinWidth()-1);
                        			lastFeature=addEnrichedDomain(currFeatures.get(cond), lastFeature, currWin, sigCounts, condCtrlCounts==null ? 0 : ctrlCounts, pval);
                        		}
                        	}
                        }
	                	currBin++;
	                }		
				}
			}
			//Trim, quantify, & properly score currFeatures before adding them to the results
			currFeatures = processDomains(currFeatures, subRegion);
			for(ExperimentCondition cond : manager.getConditions())
				results.get(cond).addAll(currFeatures.get(cond));
			return results;
		}
		
		/**
		 * Add or append a window to the list of discovered EnrichedFeatures
		 * 
		 * @param currResults
		 * @param lastFeat
		 * @param currWin
		 * @param sigWinHits
		 * @param ctrlWinHits
		 * @param score
		 * @return
		 */
		protected EnrichedFeature addEnrichedDomain(List<EnrichedFeature> currResults, EnrichedFeature lastFeat, Region currWin, float sigWinHits, float ctrlWinHits, double score){
			EnrichedFeature resFeat=null;
			
			//Is this hit close to the previously added one? If so, merge
			if(lastFeat!=null && currWin.getStrand()==lastFeat.getCoords().getStrand() && currWin.distance(lastFeat.getCoords())<=sconfig.getFeatureMergeWindow()){
				lastFeat.setCoords(strandedEventDetection ? 
						new StrandedRegion(gen, lastFeat.getCoords().getChrom(), lastFeat.getCoords().getStart(), currWin.getEnd(), lastFeat.getCoords().getStrand()) :
						new Region(gen, lastFeat.getCoords().getChrom(), lastFeat.getCoords().getStart(), currWin.getEnd()));
				if(sigWinHits>lastFeat.getSignalCount())
					lastFeat.setSignalCount(sigWinHits);
				if(ctrlWinHits>lastFeat.getControlCount())
					lastFeat.setControlCount(ctrlWinHits);
				resFeat=lastFeat;
			}else{
				EnrichedFeature feat=null;
				try {
					feat = new EnrichedFeature(currWin, null, null, sigWinHits, ctrlWinHits, score);
					currResults.add(feat);
				} catch (Exception e) {
					e.printStackTrace();
				}
				resFeat=feat;
			}return(resFeat);
		}
		
		/**
		 * ProcessDomains: 
		 *  - Trims feature coordinates back to agree with overlapping hits. 
		 *  - Counts hits in each feature, per sample 
		 * 
		 * @param currFeatures
		 * @param current region
		 * @return : Lists of EnrichedFeatures, indexed by condition
		 */
		protected Map<ExperimentCondition, List<EnrichedFeature>> processDomains(Map<ExperimentCondition,List<EnrichedFeature>> currFeatures, Region currSubRegion){
			for(ExperimentCondition currCondition : manager.getConditions()){
				for(EnrichedFeature f : currFeatures.get(currCondition)){
					Map<Sample, List<StrandedBaseCount>> fHitsPos = overlappingHits(hitsPos, f);
					Map<Sample, List<StrandedBaseCount>> fHitsNeg = overlappingHits(hitsNeg, f);
					
					//Trim the coordinates
					trimFeature(f, fHitsPos, fHitsNeg, currCondition);
					
					//Quantify the feature in each Sample and in the condition in which it was found
					quantifyFeature(f, fHitsPos, fHitsNeg, currCondition);
				}
			}
			return(currFeatures);
		}
	
		/**
		 * Trim a feature back to the coordinates of the first and last overlapping hit, if necessary.
		 * Only trims - does not allow extension past the original bounds. 
		 * @param f
		 * @param featureHits
		 * @param featureCond
		 */
		protected void trimFeature(Feature f, Map<Sample, List<StrandedBaseCount>> featureHitsPos, Map<Sample, List<StrandedBaseCount>> featureHitsNeg, ExperimentCondition featureCond){
			int minCoord=f.getCoords().getEnd();
			int maxCoord=f.getCoords().getStart();
			for(Sample s : featureCond.getSignalSamples()){
				for(int strand=0; strand<=1; strand++){
					List<StrandedBaseCount> featureHitsCurr = strand==0 ? featureHitsPos.get(s) : featureHitsNeg.get(s);
					for(StrandedBaseCount sbc : featureHitsCurr){
						int hitL = getLeft(sbc);
						int hitR = getRight(sbc);
						
						if(hitL<minCoord)
							minCoord = hitL;
						if(hitR>maxCoord)
							maxCoord = hitR;
					}
				}
			}
			minCoord = Math.max(minCoord, f.getCoords().getStart());
			maxCoord = Math.min(maxCoord+1, f.getCoords().getEnd()); //Plus one to avoid minCoord & maxCoord on the same base if there are no extensions
			f.setCoords(f.getCoords().expand(f.getCoords().getStart()-minCoord, maxCoord-f.getCoords().getEnd()));
		}
		
		/**
		 * Count the number of tags that overlap this enriched feature in every sample and in the condition in which it was discovered.
		 * In place update of counts. 
		 * @param f
		 * @param featureHitsPos
		 * @param featureHitsNeg
		 * @param featureCond
		 */
		protected void quantifyFeature(EnrichedFeature f, Map<Sample, List<StrandedBaseCount>> featureHitsPos, Map<Sample, List<StrandedBaseCount>> featureHitsNeg, ExperimentCondition featureCond){
			//Count the hits per sample
			float[] sampCountsPos = new float[manager.getSamples().size()];
			float[] sampCountsNeg = new float[manager.getSamples().size()];
			for(Sample s : manager.getSamples()){
				float currSampPos=0, currSampNeg=0;
				for(StrandedBaseCount sbc : featureHitsPos.get(s))
					currSampPos+=sbc.getCount();
				for(StrandedBaseCount sbc : featureHitsNeg.get(s))
					currSampNeg+=sbc.getCount();
				sampCountsPos[s.getIndex()]=currSampPos;
				sampCountsNeg[s.getIndex()]=currSampNeg;
			}
			
			//Sample-pooled count for the condition that this feature was discovered in
			float pooledSigCount=0, pooledCtrlCount=0;
			for(Sample s : featureCond.getSignalSamples())
				pooledSigCount+=sampCountsPos[s.getIndex()]+sampCountsNeg[s.getIndex()];
			for(Sample s : featureCond.getControlSamples())
				pooledCtrlCount+=sampCountsPos[s.getIndex()]+sampCountsNeg[s.getIndex()];
			//If there is no control, use the expected number of tags in this region as the control count for the purposes of significance testing
			float binomialTestCtrlCount = (float) (featureCond.getControlSamples().size()==0 ? 
					(featureCond.getTotalSignalCount()*f.getCoords().getWidth())/(econfig.getMappableGenomeLength())
					:(pooledCtrlCount*(float)featureCond.getPooledSampleControlScaling()));
			
			//Update the feature
			f.setSampleCountsPos(sampCountsPos);
			f.setSampleCountsNeg(sampCountsNeg);
			f.setSignalCount(pooledSigCount);
			f.setControlCount(pooledCtrlCount);
			f.setScore(stats.binomialPValue(binomialTestCtrlCount, (pooledSigCount+binomialTestCtrlCount), sconfig.getMinSigCtrlFoldDifference()));
		}
		
	}//End of DomainFinderThread

 
}//End of DomainFinder
