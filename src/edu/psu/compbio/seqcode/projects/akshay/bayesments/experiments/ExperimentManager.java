package edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments;

import java.io.File;
import java.io.IOException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.Sample;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.BEDFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.BowtieFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.ElandFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.HitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.IDXFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.NovoFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.ReadDBHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.SAMFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.TophatFileHitLoader;
import edu.psu.compbio.seqcode.projects.shaun.EventMetaMaker;

public class ExperimentManager {
	
	protected Config config;
	protected Genome gen;
	protected HashMap<String, HitLoader> loaders = new HashMap<String, HitLoader>();
	protected HashMap<String, Sample> allSamples = new HashMap<String, Sample>();
	protected List<ControlledExperiment> replicateList = new ArrayList<ControlledExperiment>();
	protected HashMap<String, ControlledExperiment> allReplicates =  new HashMap<String, ControlledExperiment>();
	protected List<ExperimentCondition> conditionList = new ArrayList<ExperimentCondition>();
	protected HashMap<String, ExperimentCondition> allConditions = new HashMap<String, ExperimentCondition>();
	protected List<ExperimentFeature> featureList = new ArrayList<ExperimentFeature>();
	protected HashMap<String, ExperimentFeature> allFeatures = new HashMap<String, ExperimentFeature>();
	protected ExperimentSet experiments;
	protected List<Point> locations;
	
	public ExperimentManager(Config conf, boolean loadReads) {
		this.config = conf;
		this.gen = this.config.getGenome();
		int repCount=0, conCount=0, samCount=0, feaCount=0;
		
		//loading the genomic locations (factor binding locations of the focus factor)
		File peaksFile = config.getPeaksFile();
		try {
			locations = EventMetaMaker.loadPoints(peaksFile, this.gen);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		// loading hitloaders in loaders hashmap
		for(ExptDescriptor e: this.config.getExperiments()){
			System.err.println("Processing HitLoaders for:\t"+e.condition+"\t"+e.replicate);
			for(Pair<String, String> src : e.sources){
				String name = src.car();
				String type = src.cdr();
				System.err.println("Loading from "+type+" hit loader:\t"+name);
				if(src.cdr().equals("READDB")){
					HitLoader hl = this.getReadDBHitLoader(name);
					loaders.put(name, hl);
				}
				else{
					HitLoader hl = getFileHitLoader(name, type, config.getNonUnique());
					loaders.put(name, hl);
				}
			}
		}
		
		//loading Samples to allSamples and sampleList
		for(ExptDescriptor e: this.config.getExperiments()){
			String sampleName;
			if(e.signal)
				sampleName = e.feature+":"+e.condition+":"+e.replicate+":signal";
			else
				sampleName = e.feature+":"+e.condition+":"+e.replicate+":control";
			if(!this.allSamples.containsKey(sampleName)){
				Sample samp = new Sample(samCount, this.config, sampleName, e.perBaseReadLimit,e.winsize);
				this.allSamples.put(sampleName, samp);
				samCount++;
			}
			for(Pair<String,String> source : e.sources){
				String name = source.car();
				this.allSamples.get(sampleName).addHitLoader(this.loaders.get(name));
				
			}
			if(loadReads)
				this.allSamples.get(sampleName).loadHits();
		}
		
		//merging the learned genome incase the the user has not provided the genome
		if(this.gen == null){
			List<Genome> estGenomes = new ArrayList<Genome>();
			for(String s : this.allSamples.keySet())
				estGenomes.add(this.allSamples.get(s).getGenome());
			this.gen = this.config.mergeGenomes(estGenomes);
			for(String s : this.allSamples.keySet())
				this.allSamples.get(s).setGenome(this.gen);
		}
		
		//initialize the replicates (scaling is estimated by default between the signal and control)
		for(ExptDescriptor e: this.config.getExperiments()){
			if(e.signal){
				String repName = e.feature+":"+e.condition+":"+e.replicate;
				if(!allReplicates.containsKey(repName)){
					Sample sig=null, ctrl=null;
					if(allSamples.containsKey(repName+":signal")){ //Require that there is a signal (in case of orphan/default controls)
						if(!loadReads)
							allSamples.get(repName+":signal").loadHits();
						sig = allSamples.get(repName+":signal");
					}
					Region[] reg= new Region[locations.size()];
					int pointCount = 0;
					for(Point p : locations){
						reg[pointCount] = new Region(this.gen,p.getChrom()
									,p.getLocation()-sig.getWinSize(),p.getLocation()+sig.getWinSize());
						pointCount++;
					}
					sig.setRegionCounts(reg);
					String cntrl_name="";	
					if(allSamples.containsKey(repName+":control")){           //Ctrl1: if there is a control defined for this condition & replicate
						if(!loadReads)
							allSamples.get(repName+":control").loadHits();
						ctrl = allSamples.get(repName+":control");
						ctrl.setRegionCounts(reg);
						cntrl_name = repName+":control";
					}
					else if(allSamples.containsKey(e.feature+":"+e.condition+":DEFAULT:control")){ //Ctrl2: if there is a default control for this condition 
						if(!loadReads)
							allSamples.get(e.feature+":"+e.condition+":DEFAULT:control").loadHits();
						ctrl = allSamples.get(e.feature+":"+e.condition+":DEFAULT:control");
						ctrl.setRegionCounts(reg);
						cntrl_name = e.feature+":"+e.condition+":DEFAULT:control";
					}
					else if(allSamples.containsKey("DEFAULT:DEFAULT:DEFAULT:control")){             //Ctrl3: if there is a global default control
						if(!loadReads)	
							allSamples.get("DEFAULT:DEFAULT:DEFAULT:control").loadHits();
						ctrl = allSamples.get("DEFAULT:DEFAULT:DEFAULT:control");
						ctrl.setRegionCounts(reg);
						cntrl_name = "DEFAULT:DEFAULT:DEFAULT:control";
					}	
					//If no control specified, ctrl is still null
						
					ControlledExperiment rep = new ControlledExperiment(this.config, repCount, e.feature, e.condition, e.replicate, 
							sig, ctrl, e.bindingmodel, this.config.getEstimateScaling(), this.config.getScalingByMedian());
					if(!loadReads){
						rep.flushReads();
						allSamples.get(repName+":signal").flushCounts();
						allSamples.get(cntrl_name).flushCounts();
					}
					
					allReplicates.put(repName, rep);
					replicateList.add(rep);
					repCount++;
				}
			}
		}
		
		
		//initializing the condition lists
		List<String> replicatesByConditionNames = new ArrayList<String>();
		List<List<ControlledExperiment>> replicatesByConditionReps = new ArrayList<List<ControlledExperiment>>();
		for(ExptDescriptor e : this.config.getExperiments()){
			String repName = e.feature+":"+e.condition+":"+e.replicate;
			if(allReplicates.containsKey(repName)){
				if(!replicatesByConditionNames.contains(e.condition)){
					replicatesByConditionReps.add(new ArrayList<ControlledExperiment>());
					replicatesByConditionNames.add(e.condition);
				}
				int index = replicatesByConditionNames.indexOf(e.condition);
				List<ControlledExperiment> currReps = replicatesByConditionReps.get(index);
				if(!currReps.contains(allReplicates.get(repName))){
					currReps.add(allReplicates.get(repName));
				}
			}
		}
		for(String s: replicatesByConditionNames){
			int index = replicatesByConditionNames.indexOf(s);
			conditionList.add(new ExperimentCondition(config, conCount, s, replicatesByConditionReps.get(index)));
			allConditions.put(s, new ExperimentCondition(config, conCount, s, replicatesByConditionReps.get(index)));
			conCount++;
		}
		
		//initializing the feature lists
		List<String> conditionsByFeatureNames = new ArrayList<String>();
		List<List<ExperimentCondition>> conditionsByFeatureCons = new ArrayList<List<ExperimentCondition>>();
		for(ExptDescriptor e : this.config.getExperiments()){
			if(!conditionsByFeatureNames.contains(e.feature)){
				conditionsByFeatureCons.add(new ArrayList<ExperimentCondition>());
				conditionsByFeatureNames.add(e.feature);
			}
			int index = conditionsByFeatureNames.indexOf(e.feature);
			List<ExperimentCondition> currCons = conditionsByFeatureCons.get(index);
			if(!currCons.contains(allConditions.get(e.condition))){
				currCons.add(allConditions.get(e.condition));
			}
		}
		for(String s: conditionsByFeatureNames){
			int index = conditionsByFeatureNames.indexOf(s);
			featureList.add(new ExperimentFeature(this.config, feaCount, s, conditionsByFeatureCons.get(index)));
			allFeatures.put(s, new ExperimentFeature(this.config, feaCount, s, conditionsByFeatureCons.get(index)));
		}
		
		//Oveall Experimentset
		this.experiments = new ExperimentSet(this.replicateList, this.conditionList, this.featureList);
	}
	
	//Accessors
	public ExperimentSet getExperimentSet(){return this.experiments;}
	public Genome getGenome(){return this.gen;}
	public List<ExperimentCondition> getChromatinConditionList(){return 
			this.allFeatures.get("CHROMATIN").getCondtionList();}
	public List<ExperimentCondition> getFacConditionList(){return 
			this.allFeatures.get("FACTOR").getCondtionList();}

	
	public HitLoader getReadDBHitLoader(String name){
		List<SeqLocator> locs = new ArrayList<SeqLocator>();
		String[] pieces = name.trim().split(";");
        if (pieces.length == 2) {
            locs.add(new SeqLocator(pieces[0], pieces[1]));
        } else if (pieces.length == 3) {
            locs.add(new SeqLocator(pieces[0], pieces[1], pieces[2]));
        } else {
            throw new RuntimeException("Couldn't parse a ChipSeqLocator from " + name);
        }
		return (new ReadDBHitLoader(gen, locs));
	}
	
	public HitLoader getFileHitLoader(String name, String format, boolean useNonUnique){
		HitLoader currReader=null;
		File file = new File(name);
		if(!file.isFile()){System.err.println("File not found: "+file.getName());System.exit(1);}
		if(format.equals("SAM") || format.equals("BAM")){
			currReader = new SAMFileHitLoader(file,useNonUnique);
		}else if(format.equals("TOPSAM")){
			currReader = new TophatFileHitLoader(file,useNonUnique);
		}else if(format.equals("ELAND")){
			currReader = new ElandFileHitLoader(file,useNonUnique);
		}else if(format.equals("NOVO")){
			currReader = new NovoFileHitLoader(file,useNonUnique);
		}else if(format.equals("BOWTIE")){
			currReader = new BowtieFileHitLoader(file,useNonUnique);
		}else if(format.equals("BED")){
			currReader = new BEDFileHitLoader(file,useNonUnique);
		}else if(format.equals("IDX")){
			currReader = new IDXFileHitLoader(file,useNonUnique);
		}else{
		    System.err.println("Unknown file format: "+format);
		    System.exit(1);
		}
		return currReader;
	}
	
	public static void main(String[] args){
		List<Integer> tert = new ArrayList<Integer>();
		tert.add(1);
		tert.add(4);
		int fet = tert.get(0);
		fet = 777;
		System.out.println(tert.get(0));
		
	}

}
