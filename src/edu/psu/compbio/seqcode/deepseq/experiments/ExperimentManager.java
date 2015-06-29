package edu.psu.compbio.seqcode.deepseq.experiments;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import edu.psu.compbio.seqcode.deepseq.hitloaders.BEDFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.BowtieFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.HitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.IDXFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.NovoFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.ReadDBHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.SAMFileHitLoader;
import edu.psu.compbio.seqcode.deepseq.hitloaders.TophatFileHitLoader;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/** 
 * ExperimentManager acts as an interface to all experiment conditions and replicates.
 * This class serves mainly to initialize the experiment tree.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentManager {

	protected ExptConfig econfig;
	protected SeqDataLoader sdloader = null;
	protected Genome gen;
	
	//Experiment tree elements
	protected HashMap<String, HitLoader> loaders = new HashMap<String,HitLoader>();
	protected List<Sample> samples = new ArrayList<Sample>();
	protected List<ControlledExperiment> replicates = new ArrayList<ControlledExperiment>();
	protected List<ExperimentCondition> conditions = new ArrayList<ExperimentCondition>();
	protected List<ExperimentTarget> targets = new ArrayList<ExperimentTarget>();
	protected List<ExperimentType> expttypes = new ArrayList<ExperimentType>();
	
	//Lookups
	protected HashMap<ExperimentCondition, Integer> conditionIndex = new HashMap<ExperimentCondition, Integer>();
	protected HashMap<Integer, ExperimentCondition> indexedCondition = new HashMap<Integer, ExperimentCondition>();
	protected HashMap<String, ExperimentCondition> namedCondition = new HashMap<String, ExperimentCondition>();

	
	/**
	 * Constructor:
	 *    Using arguments loaded by the ExptConfig, initialize (in this order):
	 *    HitLoaders, Samples, Replicates, Conditions, Targets, ExptTypes.
	 * @param c : ExptConfig
	 * @param loadReads : boolean. for some applications, reads do not have to be loaded. Use with caution. 
	 */
	public ExperimentManager(ExptConfig c){this(c, true);}
	public ExperimentManager(ExptConfig c, boolean loadReads){
		econfig = c;
		gen = econfig.getGenome();
		
		HashMap<String, Sample> allSamples = new HashMap<String, Sample>();
		HashMap<String, ControlledExperiment> allReplicates = new HashMap<String, ControlledExperiment>();
		
		List<ExptDescriptor> descriptors = econfig.getExperimentDescriptors();
		int repCount=0, condCount=0, sampCount=0, targCount=0, etypeCount=0;
		
		//Pre-step; do we need a SeqDataLoader?
		boolean makeSeqDataLoader = false;
		for(ExptDescriptor e : descriptors){
			for(Pair<String,String> source : e.sources){
				String type = source.cdr();
				if(type.equals("READDB")) 
					makeSeqDataLoader=true;
			}
		}
		if(makeSeqDataLoader)
			try {
				sdloader = new SeqDataLoader();
			} catch (SQLException e1) {
				e1.printStackTrace();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		
		//Firstly, initialize all hit loaders. 
		//This is done in a separate first pass, because it is possible (albeit unlikely)
		//that multiple conditions share the same hit loader, and you don't want to load things twice.  
		for(ExptDescriptor e : descriptors){
			if(econfig.getPrintLoadingProgress())
				System.err.println("Processing HitLoaders for:\t"+e.condition+"\t"+e.replicate);
			for(Pair<String,String> source : e.sources){
				String name = source.car();
				String type = source.cdr();
				if(type.equals("READDB")){ //ReadDB HitLoader
					HitLoader hl = makeReadDBHitLoader(sdloader, name);
					//hit loader does not have to be sourced here -- that happens in the samples part below
					loaders.put(name, hl);
				}else{  //Assume File HitLoader
					HitLoader hl = makeFileHitLoader(name, type, econfig.getNonUnique());
					//hit loader does not have to be sourced here -- that happens in the samples part below
					loaders.put(name, hl);
				}
			}
		}
		
		//Secondly, load the samples (load each sample name once)
		for(ExptDescriptor e : descriptors){
			String sampleName;
			if(e.signal)
				sampleName = e.condition+":"+e.replicate+":signal";
			else
				sampleName = e.condition+":"+e.replicate+":control";
			
			if(econfig.getPrintLoadingProgress() && loadReads)
				System.err.print("Loading data from "+sampleName);
			
			if(!allSamples.containsKey(sampleName)){
				Sample samp = new Sample(sampCount, econfig, sampleName, e.perBaseMaxReads, e.signal);
				allSamples.put(sampleName, samp);
				samples.add(samp);
				sampCount++;
			}
			for(Pair<String,String> source : e.sources){
				String name = source.car();
				allSamples.get(sampleName).addHitLoader(loaders.get(name));
			}
			if(loadReads){
				allSamples.get(sampleName).initializeCache(econfig.getCacheAllData(), econfig.getInitialCachedRegions());
				if(econfig.getPrintLoadingProgress())
					System.err.println("\tLoaded.");
			}
		}
		//Merge estimated genomes if necessary
		if(gen == null){
			List<Genome> estGenomes = new ArrayList<Genome>();
			for(String s : allSamples.keySet())
				estGenomes.add(allSamples.get(s).getGenome());
			gen = econfig.mergeEstGenomes(estGenomes);
			for(String s : allSamples.keySet())
				allSamples.get(s).setGenome(gen);
		}
		
		//Thirdly, initialize the replicates
		for(ExptDescriptor e : descriptors){
			if(e.signal){
				String repName = e.condition+":"+e.replicate;
				if(!allReplicates.containsKey(repName)){
					Sample sig=null, ctrl=null;
					if(allSamples.containsKey(repName+":signal")){ //Require that there is a signal (in case of orphan/default controls)
						sig = allSamples.get(repName+":signal");
						
						if(allSamples.containsKey(repName+":control"))           //Ctrl1: if there is a control defined for this condition & replicate
							ctrl = allSamples.get(repName+":control");
						else if(allSamples.containsKey(e.condition+":DEFAULT:control")) //Ctrl2: if there is a default control for this condition 
							ctrl = allSamples.get(e.condition+":DEFAULT:control");
						else if(allSamples.containsKey("DEFAULT:DEFAULT:control"))             //Ctrl3: if there is a global default control
							ctrl = allSamples.get("DEFAULT:DEFAULT:control");
						//If no control specified, ctrl is still null
						
						ControlledExperiment rep = new ControlledExperiment(econfig, repCount, e.condition, e.replicate, sig, ctrl);
						allReplicates.put(repName, rep);
						replicates.add(rep);
						repCount++;
					}
				}
			}
		}
		
		//Fourthly, initialize the conditions (not using Hash any more so that ordering is maintained from design file)
		List<String> replicatesByConditionNames = new ArrayList<String>();
		List<List<ControlledExperiment>> replicatesByConditionReps = new ArrayList<List<ControlledExperiment>>();
		for(ExptDescriptor e : descriptors){
			String repName = e.condition+":"+e.replicate;
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
			conditions.add(new ExperimentCondition(econfig, condCount, s, replicatesByConditionReps.get(index), econfig.getEstimateScaling()));
			condCount++;
		}
		
		//Fifthly, initialize the targets
		List<String> replicatesByTargetNames = new ArrayList<String>();
		List<List<ControlledExperiment>> replicatesByTargetReps = new ArrayList<List<ControlledExperiment>>();
		for(ExptDescriptor e : descriptors){
			String repName = e.condition+":"+e.replicate;
			if(allReplicates.containsKey(repName)){
				if(!replicatesByTargetNames.contains(e.target)){
					replicatesByTargetReps.add(new ArrayList<ControlledExperiment>());
					replicatesByTargetNames.add(e.target);
				}
				int index = replicatesByTargetNames.indexOf(e.target);
				List<ControlledExperiment> currReps = replicatesByTargetReps.get(index);
				if(!currReps.contains(allReplicates.get(repName))){
					currReps.add(allReplicates.get(repName));
				}
			}
		}
		for(String s: replicatesByTargetNames){
			int index = replicatesByTargetNames.indexOf(s);
			targets.add(new ExperimentTarget(econfig, targCount, s, replicatesByConditionReps.get(index)));
			targCount++;
		}
		
		//Sixthly, initialize the types
		List<String> replicatesByExptTypeNames = new ArrayList<String>();
		List<List<ControlledExperiment>> replicatesByExptTypeReps = new ArrayList<List<ControlledExperiment>>();
		for(ExptDescriptor e : descriptors){
			String repName = e.condition+":"+e.replicate;
			if(allReplicates.containsKey(repName)){
				if(!replicatesByExptTypeNames.contains(e.expttype)){
					replicatesByExptTypeReps.add(new ArrayList<ControlledExperiment>());
					replicatesByExptTypeNames.add(e.expttype);
				}
				int index = replicatesByExptTypeNames.indexOf(e.expttype);
				List<ControlledExperiment> currReps = replicatesByExptTypeReps.get(index);
				if(!currReps.contains(allReplicates.get(repName))){
					currReps.add(allReplicates.get(repName));
				}
			}
		}
		for(String s: replicatesByExptTypeNames){
			int index = replicatesByExptTypeNames.indexOf(s);
			expttypes.add(new ExperimentType(econfig, etypeCount, s, replicatesByConditionReps.get(index)));
			etypeCount++;
		}
		
		//Finally, index everything
		for(int i=0; i<getNumConditions(); i++){
			conditionIndex.put(conditions.get(i), i);
			indexedCondition.put(i, conditions.get(i));
			namedCondition.put(conditions.get(i).getName(), conditions.get(i));
		}
		
		if(econfig.getPrintLoadingProgress()){
			System.err.println("Loaded all experiments:");
			for(ExperimentCondition cond : getConditions()){
				System.err.println(" Condition "+cond.getName()+":\t#Replicates:\t"+cond.getReplicates().size());
				for(ControlledExperiment r : cond.getReplicates()){
					System.err.println("\tReplicate:\t"+r.getName());
					if(r.getControl()==null)
						System.err.println(String.format("\t\tSignal:\t%.1f", r.getSignal().getHitCount()));
					else
						System.err.println(String.format("\t\tSignal:\t%.1f\tControl:\t%.1f\tScalingFactor:\t%.3f",  r.getSignal().getHitCount(), r.getControl().getHitCount(), r.getControlScaling()));
				}
				if(cond.getTotalControlCount()>0)
					System.err.println(String.format("\tPooled replicates for condition:\t%s\n\t\tSignal:\t%.1f\tControl:%.1f\tScalingFactor:%.3f",cond.getName(), cond.getTotalSignalCount(), cond.getTotalControlCount(), cond.getPooledSampleControlScaling()));
				else
					System.err.println(String.format("\tPooled replicates for condition:\t%s\n\t\tSignal:\t%.1f",cond.getName(), cond.getTotalSignalCount()));
			}
		}

		if(sdloader!=null)
			sdloader.close();

	}
	
	//Accessors
	public List<Sample> getSamples(){return samples;}
	public List<ExperimentCondition> getConditions(){return conditions;}
	public List<ControlledExperiment> getReplicates(){return replicates;}
	public List<ExperimentTarget> getTargets(){return targets;}
	public List<ExperimentType> getExptTypes(){return expttypes;}
	public int getConditionIndex(ExperimentCondition c){return conditionIndex.get(c);}
	public ExperimentCondition getIndexedCondition(int index){return indexedCondition.get(index);}
	public ExperimentCondition getNamedCondition(String name){return namedCondition.get(name);}
	public int getNumConditions(){return conditions.size();}
	
	/**
	 * Add a ReadDB HitLoader.
	 * @param locs List of ChipSeqLocators
	 */
	private HitLoader makeReadDBHitLoader(SeqDataLoader loader, String name){
		if(loader==null)
			throw new RuntimeException("SeqDataLoader not initialized in makeReadDBHitLoader");
		
		List<SeqLocator> locs = new ArrayList<SeqLocator>();
		String[] pieces = name.trim().split(";");
		Set<String> reps = new TreeSet<String>(); 
		if (pieces.length == 2) {
			locs.add(new SeqLocator(pieces[0], reps, pieces[1]));
        } else if (pieces.length == 3) {
        	reps.add(pieces[1]);
        	locs.add(new SeqLocator(pieces[0], reps, pieces[2]));
        } else {
            throw new RuntimeException("Couldn't parse a SeqLocator from " + name);
        }
		return (new ReadDBHitLoader(loader, gen, locs, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs()));
	}
	

	/**
	 * Add a File HitLoader. File formats accepted include:
	 * IDX, NOVO, BOWTIE, BED	, SAM, TOPSAM
	 * @param files List of File/String Pairs, where the string is a format descriptor
	 */
	private HitLoader makeFileHitLoader(String filename, String format, boolean useNonUnique){
		HitLoader currReader=null;
		File file = new File(filename);
		if(!file.isFile()){System.err.println("File not found: "+file.getName());System.exit(1);}
		if(format.equals("SAM") || format.equals("BAM")){
			currReader = new SAMFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("TOPSAM")){
			currReader = new TophatFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("NOVO")){
			currReader = new NovoFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("BOWTIE")){
			currReader = new BowtieFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("BED")){
			currReader = new BEDFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else if(format.equals("IDX")){
			currReader = new IDXFileHitLoader(file,useNonUnique, econfig.getLoadType1Reads(), econfig.getLoadType2Reads(), econfig.getLoadPairs());
		}else{
		    System.err.println("Unknown file format: "+format);
		    System.exit(1);
		}
		return currReader;
	}
	
	

	/**
	 * Call any cleanup methods
	 */
	public void close(){
		for(String l : loaders.keySet()){
			loaders.get(l).cleanup();
		}
		for(Sample s : samples){
			s.close();
		}
	}
	
	/**
	 * This main method is only for testing the ExperimentManager system
	 * @param args
	 */
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		if(econfig.helpWanted()){
			System.err.println("ExperimentManager debugging:");
			System.err.println(econfig.getArgsList());
		}else{
			ExperimentManager manager = new ExperimentManager(econfig);
			
			System.err.println("ExptTypes:\t"+manager.getExptTypes().size());
			for(ExperimentType t : manager.getExptTypes()){
				System.err.println("ExptType "+t.getName()+":\t#Experiments:\t"+t.getExptTypeExperiments().size());
			}
			System.err.println("ExptTargets:\t"+manager.getTargets().size());
			for(ExperimentTarget t : manager.getTargets()){
				System.err.println("Target "+t.getName()+":\t#Experiments:\t"+t.getTargetExperiments().size());
			}
			
			manager.close();
		}
	}
}
