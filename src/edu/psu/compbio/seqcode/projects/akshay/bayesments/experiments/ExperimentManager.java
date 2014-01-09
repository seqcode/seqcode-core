package edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments;

import java.io.File;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.BEDFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.BowtieFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.ElandFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.HitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.IDXFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.NovoFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.ReadDBHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.SAMFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.TophatFileHitLoader;

public class ExperimentManager {
	
	protected Config config;
	protected Genome gen;
	protected HashMap<String, HitLoader> loaders = new HashMap<String, HitLoader>();
	protected List<Sample> sampleList = new ArrayList<Sample>();
	protected HashMap<String, Sample> allSamples = new HashMap<String, Sample>();
	protected List<ControlledExperiment> replicateList = new ArrayList<ControlledExperiment>();
	protected HashMap<String, ControlledExperiment> allReplicates =  new HashMap<String, ControlledExperiment>();
	protected List<ExperimentCondition> allConditions = new ArrayList<ExperimentCondition>();
	protected HashMap<String, ExperimentCondition> conditionList = new HashMap<String, ExperimentCondition>();
	protected List<ExperimentFeature> allFeatures = new ArrayList<ExperimentFeature>();
	protected HashMap<String, ExperimentFeature> featureList = new HashMap<String, ExperimentFeature>();
	protected ExperimentSet experiments;
	
	public ExperimentManager(Config conf, boolean loadReads) {
		this.config = conf;
		this.gen = this.config.getGenome();
		int repCount=0, conCount=0, samCount=0, feaCount=0;
		
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
			
		}
		
		
		
	}
	
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
	

}
