package org.seqcode.projects.seed;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.projects.seed.features.EnrichedFeature;
import org.seqcode.projects.seed.features.Feature;
import org.seqcode.projects.seed.features.SuperEnrichedFeature;



/**
 * 
 * @author akshaykakumanu
 *
 */
public class SuperEnhancerFinder extends DomainFinder{
	
	// pseudocount to calculate fold-change while finding inflectioin point for super enhancers
	public final int PSEUDOCOUINT  = 1;
	
	// At the moment only one TSS is alowed to overlap identified super enahncers
	public final int NO_OVERLAPPING_TSS_ALLOWED = 1;
	
	// P-value threshold for selcting enriched features to be stitched 
	public final double PER_ENRICHED_FEATURE_THRES = 0.001;
	
	// Minimum length of an enriched feature to be considered for stitching to get a super enhancer
	public final int MIN_ENRICHED_FEATURE_LENGHT = 100;
	
	//Srep size to calcluate slope
	public final int SLOPE_CALCULATING_STEP_SIZE = 40;
	
	//Min no of data points to change at inflection point
	public final double GRADIENT_STRICT = 1.0;
	public final double GRADIENT_RELAX = 0.9;
	
	
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
		return "org.seqcode.projects.seed.SuperEnhancerFinder";
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
	
	
	public Map<ExperimentCondition, List<Feature>> postProcess() {
		System.err.println("\nSuper Enhancer finding complete.");
		
		
		//Sort super enhancers by signal count
		
		Map<ExperimentCondition, List<SuperEnrichedFeature>> signal_sorted_features = new HashMap<ExperimentCondition,List<SuperEnrichedFeature>>();
		for( ExperimentCondition ec : manager.getConditions()){
			signal_sorted_features.put(ec, new ArrayList<SuperEnrichedFeature>());
			for(Feature f : features.get(ec)){
				SuperEnrichedFeature sfe = (SuperEnrichedFeature)f;
				sfe.setScore(sfe.getSignalCount());
				signal_sorted_features.get(ec).add(sfe);
			}
			Collections.sort(signal_sorted_features.get(ec));
		}
		
		//Calculating slope to find inflection point
		int inflection_point_index_strict = 0;
		int inflection_point_index_relax = 0;
		
		for( ExperimentCondition ec : manager.getConditions()){
			boolean reachedInflectionPoint_strict = false;
			boolean reachedInflectionPoint_relax = false;
			for(int i=0; i<signal_sorted_features.get(ec).size(); i++){
				int x=i-(int)(SLOPE_CALCULATING_STEP_SIZE/2);
				int y = i+(int)(SLOPE_CALCULATING_STEP_SIZE/2);
				if(x<0){x=0;}
				if(y >signal_sorted_features.get(ec).size()){y=signal_sorted_features.get(ec).size();}
				
			
				float[] vals = new float[y-x];
				for(int j=x;j<y;j++){
					vals[j-x] = signal_sorted_features.get(ec).get(j).getSignalCount();
				}
				
				double slope = getSlope(vals);
				
				if(slope >=GRADIENT_STRICT && !reachedInflectionPoint_strict){
					inflection_point_index_strict = i;
					reachedInflectionPoint_strict = true;
				}
				
				if(slope >=GRADIENT_RELAX && !reachedInflectionPoint_relax){
					inflection_point_index_relax = i;
					reachedInflectionPoint_relax = true;
				}
				
				signal_sorted_features.get(ec).get(i).setSlope(slope);
			}
		}
		
		Map<ExperimentCondition, List<Feature>> finalfeature = new HashMap<ExperimentCondition,List<Feature>>();
		
		for(ExperimentCondition ec : manager.getConditions()){
			finalfeature.put(ec, new ArrayList<Feature>());
			finalfeature.get(ec).addAll(signal_sorted_features.get(ec));
		}
		
		Map<ExperimentCondition, List<Feature>> superEnhancers_strict = new HashMap<ExperimentCondition,List<Feature>>();
		Map<ExperimentCondition, List<Feature>> superEnhancers_relax = new HashMap<ExperimentCondition,List<Feature>>();
		
		for(ExperimentCondition ec : manager.getConditions()){
			superEnhancers_strict.put(ec, new ArrayList<Feature>());
			superEnhancers_relax.put(ec, new ArrayList<Feature>());
			int index = 0;
			for(Feature ff :finalfeature.get(ec)){
				if(index >=  inflection_point_index_strict){
					superEnhancers_strict.get(ec).add(ff);
				}
				if(index >= inflection_point_index_relax){
					superEnhancers_relax.get(ec).add(ff);
				}
				index++;
			}
			Collections.reverse(superEnhancers_strict.get(ec));
			Collections.reverse(superEnhancers_relax.get(ec));
		}
	
		
		//All features
		this.printEventsFile(finalfeature, ".all.domains");

		// Super enhancer features
		this.printEventsFile(superEnhancers_strict, ".superenhancer_strict.domains");
		this.printEventsFile(superEnhancers_relax, ".superenhancer_relax.domains");
		
		//TFs of SFs
		this.printTFsOfSFs(superEnhancers_relax, ".TFs_of_SFs_relax.domains");
		this.printTFsOfSFs(superEnhancers_strict, ".TFs_of_SFs_strict.domains");
		
		return features;
	}
	
	/**
	 * Help String for Super Enhancer Finder
	 * @return
	 */
	public static String getDomainFinderArgs() {
		return(new String("" +
				"SuperEnhancerFinder arguments:\n"+
				"\t--supStitchWin <Window to stitch typical enhancers into super enhancers (default: supStitchWin=12500bp)>\n" +
				"\t--distalDistance <Minmum distnace from TSSs for calling domains (default: distalDistance= 2000bp)>\n" +
				"\t--refTSSs <TSS annotations, eg :- chr2:45667-45767:+ >\n" +
				""));
	}
	
	protected double getSlope(float values[]){
		double slope=0;
		for(int i=0; i<(values.length-1); i++){
			if(values[i+1] >values[i]){
				slope++;
			}
		}
		return slope/(values.length-1);
	}
	
	protected void printTFsOfSFs(Map<ExperimentCondition,List<Feature>> SuperFeatures,String suffix){
		try {
			for(ExperimentCondition cond : manager.getConditions()){
				String condName = cond.getName().equals("experiment") ? "" : "_"+cond.getName();
	    		String filename = sconfig.getOutputParentDir()+File.separator+sconfig.getOutBase()+condName+suffix;
				FileWriter fout = new FileWriter(filename);
				fout.write("#"+getProgramName()+"\n");
				fout.write("#Arguments:\t"+sconfig.getArgs()+"\n");
				fout.write("##date "+(new Date()).toString()+"\n");
				for(Feature f : SuperFeatures.get(cond)){
					if(f instanceof SuperEnrichedFeature){
						SuperEnrichedFeature sf = (SuperEnrichedFeature) f;
						for(EnrichedFeature ef : sf.getTEFs()){
							fout.write(ef.getCoords().getLocationString()+"\t"+sf.getCoords().getLocationString()+"\n");
						}
					}
				}
				fout.close();
			}
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("SuperEnhancerFinder version "+SuperEnhancerFinder.version+"\n\n");
		
		String SuperEnhancerFinderHelp = new String("" +
				"SuperEnhancerFinder arguments:\n"+
				"\t--supStitchWin <Window to stitch typical enhancers into super enhancers (default: supStitchWin=12500bp)>\n" +
				"\t--distalDistance <Minmum distnace from TSSs for calling domains (default: distalDistance= 2000bp)>\n" +
				"\t--refTSSs <TSS annotations, eg :- chr2:45667-45767:+ >\n" +
				"");
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		SEEDConfig scon = new SEEDConfig(gcon, args);
		
		if(scon.helpWanted()){
			//System.out.println(DomainFinder.getDomainFinderArgs());
			System.err.println(gcon.getArgsList()+
					econ.getArgsList()+
					scon.getArgsList()+
					SuperEnhancerFinderHelp);
		}else{
			ExperimentManager man = new ExperimentManager(econ);
			SuperEnhancerFinder finder = new SuperEnhancerFinder(gcon, econ, scon, man);
			System.err.println("\nBeginning Super Enhancer finding...");
			finder.execute();
			man.close();
			
		}
	}
	
	
	
	/**
	 * SuperEnhancetFindetThread that searches for super enhancers
	 * @author akshaykakumanu
	 *
	 */
	public class SuperEnhancerFinderThread extends DomainFinderThread{

		public SuperEnhancerFinderThread(List<Region> regs) {
			super(regs);
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
			
			
			// Remove enriched domains less than 100bp. Playing around with the parameter at the moment (Needed especially when extending reads)
			// Also removes enriched feature that are have a p-value < 0.001(PER_ENRICHED_FEATURE_THRES). Playing around with the parameter at the moment
			for (ExperimentCondition ec : manager.getConditions()){
				Iterator<EnrichedFeature> it = currFeatures.get(ec).iterator();
				while(it.hasNext()){
					EnrichedFeature ef = it.next();
					Map<Sample, List<StrandedBaseCount>> efHitsPos = overlappingHits(hitsPos, ef);
					Map<Sample, List<StrandedBaseCount>> efHitsNeg = overlappingHits(hitsNeg, ef);
					quantifyFeature(ef,efHitsPos,efHitsNeg,ec);
					if(ef.getCoords().getWidth() < MIN_ENRICHED_FEATURE_LENGHT || ef.getScore() >= PER_ENRICHED_FEATURE_THRES){
						it.remove();
					}
				}
			}
			
			
			// Do a bunch of more operations of the enriched features
			for(ExperimentCondition ec :  manager.getConditions()){
				// Stitching typical enhancers into super enhancers
				List<SuperEnrichedFeature> stitchedTPEs = this.stitchTypicalEnhancers(currFeatures.get(ec));
				
				// Removing SEs spanning more than 1 TSSs
				removeMultiGeneSpanningSEs(stitchedTPEs);
				
				// Quantify the features 
				for(SuperEnrichedFeature sfe : stitchedTPEs){
					Map<Sample, List<StrandedBaseCount>> sfeHitsPos = overlappingHits(hitsPos, sfe);
					Map<Sample, List<StrandedBaseCount>> sfeHitsNeg = overlappingHits(hitsNeg, sfe);
					
					//Trim the coordinates
					trimFeature(sfe, sfeHitsPos, sfeHitsNeg, ec);
					
					//Quantify the feature in each Sample and in the condition in which it was found
					quantifyFeature(sfe,sfeHitsPos,sfeHitsNeg,ec);
				}
				
				superEnhancers.get(ec).addAll(stitchedTPEs);
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
						int dis = Math.min(sp.distance(ef.getCoords().startPoint()),sp.distance(ef.getCoords().endPoint()));
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
		
		protected List<SuperEnrichedFeature> stitchTypicalEnhancers(List<EnrichedFeature> enhancers){
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
				if(numTSSs > NO_OVERLAPPING_TSS_ALLOWED){
					it.remove();
				}
			}
			
		}
			
		
	}

}
