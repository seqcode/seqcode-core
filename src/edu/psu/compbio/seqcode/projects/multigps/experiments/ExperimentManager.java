package edu.psu.compbio.seqcode.projects.multigps.experiments;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.BEDFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.BowtieFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.ElandFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.HitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.IDXFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.NovoFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.ReadDBHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.SAMFileHitLoader;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.TophatFileHitLoader;

/** 
 * ExperimentManager acts as an interface to all experiment conditions and replicates.
 * This class serves mainly to initialize the experiment tree.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentManager {

	protected Config config;
	protected Genome gen;
	protected HashMap<String, HitLoader> loaders = new HashMap<String,HitLoader>();
	protected HashMap<String, Sample> allSamples = new HashMap<String, Sample>();
	protected HashMap<String, ControlledExperiment> allReplicates = new HashMap<String, ControlledExperiment>();
	protected List<ControlledExperiment> replicateList = new ArrayList<ControlledExperiment>();
	protected List<ExperimentCondition> conditionList = new ArrayList<ExperimentCondition>();
	protected ExperimentSet experiments = null;
	protected List<BindingEvent> events = new ArrayList<BindingEvent>();
	
	
	/**
	 * Constructor:
	 *    Using arguments loaded by the ArgsHandler, initialize (in this order):
	 *    HitLoaders, Samples, Replicates, Conditions, The ExperimentSet.
	 * @param c : Config
	 * @param loadReads : boolean. for some applications, reads do not have to be loaded. Use with caution. 
	 */
	public ExperimentManager(Config c){this(c, true);}
	public ExperimentManager(Config c, boolean loadReads){
		config = c;
		gen = config.getGenome();
		List<ExptDescriptor> descriptors = config.getExperiments();
		int repCount=0, condCount=0, sampCount=0;
		
		//Firstly, initialize all hit loaders. 
		//This is done in a separate first pass, because it is possible (albeit unlikely)
		//that multiple conditions share the same hit loader, and you don't want to load things twice.  
		for(ExptDescriptor e : descriptors){
			System.err.println("Processing HitLoaders for:\t"+e.condition+"\t"+e.replicate);
			for(Pair<String,String> source : e.sources){
				String name = source.car();
				String type = source.cdr();
				System.err.println("Loading from "+type+" hit loader:\t"+name);
				if(type.equals("READDB")){ //ReadDB HitLoader
					HitLoader hl = getReadDBHitLoader(name);
					//hit loader does not have to be sourced here -- that happens in the samples part below
					loaders.put(name, hl);
				}else{  //Assume File HitLoader
					HitLoader hl = getFileHitLoader(name, type, config.getNonUnique());
					//hit loader does not have to be sourced here -- that happens in the samples part below
					loaders.put(name, hl);
				}
			}
		}
		
		//Secondly, load the samples
		for(ExptDescriptor e : descriptors){
			String sampleName;
			if(e.signal)
				sampleName = e.condition+":"+e.replicate+":signal";
			else
				sampleName = e.condition+":"+e.replicate+":control";
			if(!allSamples.containsKey(sampleName)){
				Sample samp = new Sample(sampCount, config, sampleName, e.perBaseMaxReads);
				allSamples.put(sampleName, samp);
				sampCount++;
			}
			for(Pair<String,String> source : e.sources){
				String name = source.car();
				allSamples.get(sampleName).addHitLoader(loaders.get(name));
			}
			if(loadReads)
				allSamples.get(sampleName).loadHits();
		}
		//Merge estimated genomes if necessary
		if(gen == null){
			List<Genome> estGenomes = new ArrayList<Genome>();
			for(String s : allSamples.keySet())
				estGenomes.add(allSamples.get(s).getGenome());
			gen = config.mergeGenomes(estGenomes);
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
						
						ControlledExperiment rep = new ControlledExperiment(config, repCount, e.condition, e.replicate, sig, ctrl, e.bindingModel, config.getEstimateScaling());
						allReplicates.put(repName, rep);
						replicateList.add(rep);
						repCount++;
					}
				}
			}
		}
		
		//Fourthly, initialize the conditions (not using Hash any more so that ordering is maintained from design file)
		//HashMap<String, List<ControlledExperiment>> replicatesByCondition=new HashMap<String, List<ControlledExperiment>>();
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
			conditionList.add(new ExperimentCondition(config, condCount, s, replicatesByConditionReps.get(index)));
			condCount++;
		}
		
		
		//Finally, the overall ExperimentSet
		experiments = new ExperimentSet(conditionList, replicateList);
	}
	
	//Accessors
	public ExperimentSet getExperimentSet(){return experiments;}
	public int getNumConditions(){return experiments.getConditions().size();}
	public List<BindingEvent> getEvents(){return events;}
	public void setEvents(List<BindingEvent> e){events =e;}
	
	/**
	 * Get the maximum model width
	 * @return
	 */
	public int getMaxModelRange(){
		int max=0;
		for(String s : allReplicates.keySet())
			if(allReplicates.get(s).getBindingModel().getInfluenceRange()>max)
				max = allReplicates.get(s).getBindingModel().getInfluenceRange();
		return max;
	}
	
	/**
	 * Add a ReadDB HitLoader.
	 * @param locs List of ChipSeqLocators
	 */
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
	

	/**
	 * Add a File HitLoader. File formats accepted include:
	 * ELAND, NOVO, BOWTIE, BED	, SAM, TOPSAM
	 * @param files List of File/String Pairs, where the string is a format descriptor
	 */
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
	
	/**
	 * For each controlled experiment, simply calculate the proportion of reads in the provided 
	 * list of binding events to everything else. 
	 * @param regs
	 */
	public void estimateSignalProportion(List<BindingEvent> signalEvents){
		for(ExperimentCondition c : conditionList){
			for(ControlledExperiment r : c.getReplicates()){
				double repSigCount =0, repNoiseCount=0;
				for(BindingEvent event : signalEvents){
					if(event.isFoundInCondition(c))
						repSigCount += event.getRepSigHits(r);
				}
				repNoiseCount = r.getSignal().getHitCount() - repSigCount;
				r.setSigNoiseCounts(repSigCount,  repNoiseCount);
				System.err.println(r.getName()+"\t"+r.getIndex()+"\tsignal-noise ratio:\t"+String.format("%.4f",r.getSigProp()));
			}
		}
	}
	/**
	 * Count the binding events present in a given condition
	 * @param cond
	 * @return
	 */
	public int countEventsInCondition(ExperimentCondition cond){
		int count=0;
		for(BindingEvent e : events){
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond) <=config.getQMinThres())
				count++;
		}
		return count;
	}
	/**
	 * Count the differential binding events present in a given pair of conditions
	 * @param cond
	 * @return
	 */
	public int countDiffEventsBetweenConditions(ExperimentCondition cond, ExperimentCondition othercond){
		int count=0;
		for(BindingEvent e : events){
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond) <=config.getQMinThres())
    			if(e.getInterCondP(cond, othercond)<=config.getDiffPMinThres() && e.getInterCondFold(cond, othercond)>0)
    				count++;
		}
		return count;
	}
    /**
     * Print all binding events to files
     */
    public void writeBindingEventFiles(String filePrefix){
    	if(events.size()>0){
    		
	    	try {
	    		//Full output table (all non-zero components)
	    		String filename = filePrefix+".all.events.table";
	    		FileWriter fout = new FileWriter(filename);
	    		fout.write(BindingEvent.fullHeadString()+"\n");
	    		for(BindingEvent e : events)
	    			fout.write(e.toString()+"\n");
				fout.close();
	    		
	    		//Per-condition event files
	    		for(ExperimentCondition cond : experiments.getConditions()){
	    			//Sort on the current condition
	    			BindingEvent.setSortingCond(cond);
	    			Collections.sort(events, new Comparator<BindingEvent>(){
	    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
	    	        });
	    			//Print events
	    			String condName = cond.getName(); 
	    			condName = condName.replaceAll("/", "-");
	    			filename = filePrefix+"_"+condName+".events";
					fout = new FileWriter(filename);
					fout.write(BindingEvent.conditionHeadString(cond)+"\n");
			    	for(BindingEvent e : events){
			    		double Q = e.getCondSigVCtrlP(cond);
			    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
			    		if(e.isFoundInCondition(cond) && Q <=config.getQMinThres())
			    			fout.write(e.getConditionString(cond)+"\n");
			    	}
					fout.close();
	    		}
	    		
	    		//Differential event files
	    		if(getNumConditions()>1 && config.getRunDiffTests()){
	    			for(ExperimentCondition cond : experiments.getConditions()){
		    			//Sort on the current condition
		    			BindingEvent.setSortingCond(cond);
		    			Collections.sort(events, new Comparator<BindingEvent>(){
		    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
		    	        });
		    			
		    			for(ExperimentCondition othercond : experiments.getConditions()){
		    				if(!cond.equals(othercond)){
				    			//Print diff events
				    			String condName = cond.getName(); 
				    			String othercondName = othercond.getName(); 
				    			condName = condName.replaceAll("/", "-");
				    			filename = filePrefix+"_"+condName+"_gt_"+othercondName+".diff.events";
								fout = new FileWriter(filename);
								fout.write(BindingEvent.conditionShortHeadString(cond)+"\n");
						    	for(BindingEvent e : events){
						    		double Q = e.getCondSigVCtrlP(cond);
						    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
						    		if(e.isFoundInCondition(cond) && Q <=config.getQMinThres()){
						    			if(e.getInterCondP(cond, othercond)<=config.getDiffPMinThres() && e.getInterCondFold(cond, othercond)>0){
						    				fout.write(e.getConditionString(cond)+"\n");
						    			}
						    		}
						    	}
								fout.close();
		    				}
		    			}
		    		}
	    		}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
 
    /**
     * Print all binding events to files
     */
    public void writeFullEventFile(String filename){
    	if(events.size()>0){
    		try {
	    		//Full dataset table
		    	FileWriter fout = new FileWriter(filename);
		    	fout.write(BindingEvent.fullHeadString()+"\n");
		    	for(BindingEvent e : events)
		    		fout.write(e.toString()+"\n");
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
    
    /**
     * Print all motifs to file
     */
    public void writeMotifFile(String filename){
    	if(config.getFindingMotifs()){
    		try {
	    		//Full dataset table
		    	FileWriter fout = new FileWriter(filename);
		    	for(ExperimentCondition cond : experiments.getConditions()){
		    		if(cond.getFreqMatrix()!=null)
		    			fout.write(WeightMatrix.printTransfacMatrix(cond.getFreqMatrix(), cond.getName()));
		    	}
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
 
    
    /**
     * Print replicate data counts to a file
     */
    public void writeReplicateCounts(String filename){
    	if(events.size()>0){
    		try {
	    		//Full dataset table
	    		FileWriter fout = new FileWriter(filename);
	    		fout.write(BindingEvent.repCountHeadString()+"\n");
	    		for(BindingEvent e : events)
	    			fout.write(e.getRepCountString()+"\n");
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
    /**
     * Print all binding events to screen
     * TESTING ONLY
     */
    public void printBindingEvents(){
    	System.err.println(events.size()+" events found");
    	for(BindingEvent e : events){
    		System.err.println(e.toString());
    	}
    }

	/**
	 * Call any cleanup methods
	 */
	public void close(){
		for(String l : loaders.keySet())
			loaders.get(l).cleanup();
	}
	
	/**
	 * This main method is only for testing the ExperimentManager system
	 * @param args
	 */
	public static void main(String[] args){
		
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("ExperimentManager debugging:");
			System.err.println(config.getArgsList());
		}else{
			ExperimentManager manager = new ExperimentManager(config);
			
			ExperimentSet eset = manager.getExperimentSet();
			System.err.println("Conditions:\t"+eset.getConditions().size());
			for(ExperimentCondition c : eset.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : eset.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			manager.close();
		}
	}
}
