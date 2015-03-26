package edu.psu.compbio.seqcode.projects.seed;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;


import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.projects.seed.FeatureDetection.FeatureDetectionThread;
import edu.psu.compbio.seqcode.projects.seed.features.EnrichedFeature;
import edu.psu.compbio.seqcode.projects.seed.features.Feature;
import edu.psu.compbio.seqcode.projects.seed.features.SuperEnrichedFeature;

/**
 * 
 * @author akshaykakumanu
 *
 */
public class SuperEnhancerFinder extends DomainFinder{
	
	// Min distance from known TSSs for a feature to be called distal
	protected int minDistalDistance;
	
	// Super enhancer stitch distance
	protected int superStitchWin;
	
	// Reference TSSs loaded from strandedpoints peaks file
	protected List<StrandedPoint> refTSSs;
	
	public static String version = "0.1";
	
	/**
	 * Return the class name
	 */
	public String getProgramName(){
		return "edu.psu.compbio.seqcode.projects.seed.SuperEnhancerFinder";
	}
	
	/**
	 * Return a thread for this implementation
	 */
	public FeatureDetectionThread getMyThread(List<Region> regs){
		return new SuperEnhancerFinderThread(regs);
	}

	public SuperEnhancerFinder(GenomeConfig gcon, ExptConfig econ,
			SEEDConfig scon, ExperimentManager man) {
		super(gcon, econ, scon, man);
		// TODO Auto-generated constructor stub
		
		minDistalDistance = sconfig.getMinDistalDistance();
		superStitchWin = sconfig.getSuperStitchWin();
		refTSSs = sconfig.getRefTSSs();
	}
	
	
	/**
	 * SuperEnhancetFindetThread that searches for super enhancers
	 * @author akshaykakumanu
	 *
	 */
	public class SuperEnhancerFinderThread extends DomainFinderThread{

		public SuperEnhancerFinderThread(List<Region> regs) {
			super(regs);
			// TODO Auto-generated constructor stub
		}
		
		
		/**
		 * findFeatures: find potentially enriched domains on the genome by comparing tag counts to a background model,
		 * and then finds super enhancers as described in ...et.al
		 */
		public Map<ExperimentCondition, List<Feature>> findFeatures(Region subRegion){
			// Get the enriched domains using the functionality from the superclass (DomainFinderThread)
			// and the over-ridden processDomain methods below 
			return super.findFeatures(subRegion);
		}
		
		
		/**
		 * 
		 */
		protected Map<ExperimentCondition, List<EnrichedFeature>> processDomains(Map<ExperimentCondition,List<EnrichedFeature>> currFeatures, Region currSubRegion){
			Map<ExperimentCondition, List<EnrichedFeature>> superEnhancers = new HashMap<ExperimentCondition, List<EnrichedFeature>>();
			for (ExperimentCondition ec : manager.getConditions()){
				superEnhancers.put(ec, new ArrayList<EnrichedFeature>());
			}
			
			// Removing all enriched domains close to promoters
			for (ExperimentCondition ec : manager.getConditions()){
				keepOnlyDistalFeature(currFeatures.get(ec));
			}
			
			
			for(ExperimentCondition ec :  manager.getConditions()){
				// Stitching typical enhancers into super enhancers
				List<SuperEnrichedFeature> stitchedTPEs = this.stichTypicalEnhancers(currFeatures.get(ec));
				
				// Removing SEs spanning more than 2 TSSs
				removeMultiGeneSpanningSEs(stitchedTPEs);
				
				//Trim the coordinates
				for(SuperEnrichedFeature sfe: stitchedTPEs){
					trimFeature(sfe, hitsPos, hitsNeg, ec);
				}
				
				// Quantify the features 
				for(SuperEnrichedFeature sfe : stitchedTPEs){
					quantifyFeature(sfe,hitsPos,hitsNeg,ec);
				}
				
				
				
				
				
			}
			
			
			
			
			return superEnhancers;
		}
		
		/*
		 * Keep only those enriched domains that are distal to known TSSs
		 */
		protected void keepOnlyDistalFeature(List<EnrichedFeature> enrichedFeatures){
			Map<String, List<StrandedPoint>> refTSSsbyChrom = new HashMap<String, List<StrandedPoint>>();
			for(StrandedPoint sp: refTSSs){
				if(refTSSsbyChrom.containsKey(sp.getChrom())){refTSSsbyChrom.get(sp.getChrom()).add(sp);}
				else{refTSSsbyChrom.put(sp.getChrom(),new ArrayList<StrandedPoint>()); refTSSsbyChrom.get(sp.getChrom()).add(sp);}
			}
			Iterator<EnrichedFeature> it = enrichedFeatures.iterator();
			while(it.hasNext()){
				EnrichedFeature ef = it.next();
				int minDistance = Integer.MAX_VALUE;
				if(refTSSsbyChrom.containsKey(ef.getCoords().getChrom())){
					for(StrandedPoint sp : refTSSsbyChrom.get(ef.getCoords().getChrom())){
						int dis = sp.distance(ef.getCoords().getMidpoint());
						if(dis < minDistance){
							minDistance = dis;
						}
					}
				}
				if(minDistance <= minDistalDistance ){
					it.remove();
				}
			}
		}
		
		protected List<SuperEnrichedFeature> stichTypicalEnhancers(List<EnrichedFeature> enhancers){
			List<SuperEnrichedFeature> superEnhancers = new ArrayList<SuperEnrichedFeature>();
			SuperEnrichedFeature lastadded = null;
			try {
				for(EnrichedFeature ef : enhancers){
					if(lastadded == null){
						superEnhancers.add(new SuperEnrichedFeature(ef));
						lastadded = superEnhancers.get(superEnhancers.size()-1);
					}else{
						if(lastadded.getLastTEF().getCoords().getMidpoint().distance(ef.getCoords().getMidpoint()) < superStitchWin){
							superEnhancers.get(superEnhancers.size()-1).addTypicalFeature(ef);
							lastadded = superEnhancers.get(superEnhancers.size()-1);
						}else{
							superEnhancers.add(new SuperEnrichedFeature(ef));
							lastadded = superEnhancers.get(superEnhancers.size()-1);
						}
					}
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return superEnhancers;
		}
		
		protected void removeMultiGeneSpanningSEs(List<SuperEnrichedFeature> SEs){
			Map<String, List<StrandedPoint>> refTSSsbyChrom = new HashMap<String, List<StrandedPoint>>();
			for(StrandedPoint sp: refTSSs){
				if(refTSSsbyChrom.containsKey(sp.getChrom())){refTSSsbyChrom.get(sp.getChrom()).add(sp);}
				else{refTSSsbyChrom.put(sp.getChrom(),new ArrayList<StrandedPoint>()); refTSSsbyChrom.get(sp.getChrom()).add(sp);}
			}
			Iterator<SuperEnrichedFeature> it = SEs.iterator();
			while(it.hasNext()){
				int numTSSs = 0;
				SuperEnrichedFeature sef = it.next();
				if(refTSSsbyChrom.containsKey(sef.getCoords().getChrom())){
					for(StrandedPoint sp :  refTSSs){
						if(sef.getCoords().contains(sp)){
							numTSSs++;
						}
					}
				}
				if(numTSSs >=2){
					it.remove();
				}
			}
			
		}
		
		
		
		
		
	}

}
