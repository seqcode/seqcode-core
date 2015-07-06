package edu.psu.compbio.seqcode.projects.akshay.bayesments.framework;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.hitloaders.BEDFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.BowtieFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.HitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.IDXFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.NovoFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.ReadDBHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.SAMFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.TophatFileHitLoader;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.EventMetaMaker;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.BayesmentsExptDescriptor;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentFeature;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;


public class BayesmentsEMan {
	protected ExperimentManager manager;
	protected BayesmentsConfig config;
	protected Genome gen;
	protected HashMap<String, ControlledExperiment> allReplicates =  new HashMap<String, ControlledExperiment>();
	protected List<ExperimentCondition> conditionList = new ArrayList<ExperimentCondition>();
	protected HashMap<String, ExperimentCondition> allConditions = new HashMap<String, ExperimentCondition>();
	protected List<ExperimentFeature> featureList = new ArrayList<ExperimentFeature>();
	protected HashMap<String, ExperimentFeature> allFeatures = new HashMap<String, ExperimentFeature>();
	protected List<Point> locations;
	protected Map<String, Integer> condSizes = new HashMap<String, Integer>();
	
	public BayesmentsEMan(BayesmentsConfig conf, ExperimentManager em) {
		manager = em;
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
		
		
		//initializing the condition lists
		List<String> replicatesByConditionNames = new ArrayList<String>();
		List<List<ControlledExperiment>> replicatesByConditionReps = new ArrayList<List<ControlledExperiment>>();
		for(BayesmentsExptDescriptor e : this.config.getExperiments()){
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
		
		
		//initializing the feature lists
		List<String> conditionsByFeatureNames = new ArrayList<String>();
		List<List<ExperimentCondition>> conditionsByFeatureCons = new ArrayList<List<ExperimentCondition>>();
		for(BayesmentsExptDescriptor e : this.config.getExperiments()){
			if(!conditionsByFeatureNames.contains(e.feature)){
				conditionsByFeatureCons.add(new ArrayList<ExperimentCondition>());
				conditionsByFeatureNames.add(e.feature);
			}
			int index = conditionsByFeatureNames.indexOf(e.feature);
			List<ExperimentCondition> currCons = conditionsByFeatureCons.get(index);
			if(!currCons.contains(allConditions.get(e.condition))){
				currCons.add(allConditions.get(e.condition));
			}
			condSizes.put(e.condition, e.winsize);
		}
		for(String s: conditionsByFeatureNames){
			int index = conditionsByFeatureNames.indexOf(s);
			featureList.add(new ExperimentFeature(this.config, feaCount, s, conditionsByFeatureCons.get(index)));
			allFeatures.put(s, new ExperimentFeature(this.config, feaCount, s, conditionsByFeatureCons.get(index)));
		}
		
	}
	
	//Accessors
	public Genome getGenome(){return this.gen;}
	public List<ExperimentCondition> getChromatinConditionList(){return 
			this.allFeatures.get("CHROMATIN").getCondtionList();}
	public List<ExperimentCondition> getFacConditionList(){return 
			this.allFeatures.get("FACTOR").getCondtionList();}
	public Integer getConditionWinSize(String condition){return condSizes.get(condition);}
	
	
	public static void main(String[] args){
		List<Integer> tert = new ArrayList<Integer>();
		tert.add(1);
		tert.add(4);
		int fet = tert.get(0);
		fet = 777;
		System.out.println(tert.get(0));
		
	}

}
