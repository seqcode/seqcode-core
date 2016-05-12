package edu.psu.compbio.seqcode.projects.multigps;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModelPerBase;
import edu.psu.compbio.seqcode.projects.multigps.framework.MultiGPSConfig;
import edu.psu.compbio.seqcode.projects.multigps.framework.EnrichmentSignificance;
import edu.psu.compbio.seqcode.projects.multigps.framework.OutputFormatter;
import edu.psu.compbio.seqcode.projects.multigps.framework.PotentialRegionFilter;
import edu.psu.compbio.seqcode.projects.multigps.mixturemodel.BindingMixture;
import edu.psu.compbio.seqcode.projects.multigps.stats.CountsDataset;
import edu.psu.compbio.seqcode.projects.multigps.stats.DifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.EdgeRDifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.Normalization;
import edu.psu.compbio.seqcode.projects.multigps.stats.TMMNormalization;
import edu.psu.compbio.seqcode.projects.multigps.utilities.EventsPostAnalysis;

public class MultiGPS {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected MultiGPSConfig mgpsconfig;
	protected BindingManager bindingManager;
	protected BindingMixture mixtureModel;
	protected PotentialRegionFilter potentialFilter;
	protected OutputFormatter outFormatter;
	private CountsDataset data;
	protected Normalization normalizer;
	protected Map<ControlledExperiment, List<BindingModel>> repBindingModels;
	
	public MultiGPS(GenomeConfig gcon, ExptConfig econ, MultiGPSConfig c, ExperimentManager eMan){
		gconfig = gcon;
		econfig = econ;
		manager = eMan;
		mgpsconfig = c;
		mgpsconfig.makeGPSOutputDirs(true);
		outFormatter = new OutputFormatter(mgpsconfig);
		bindingManager = new BindingManager(mgpsconfig, manager);
		
		//Initialize binding models & binding model record
		repBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();
		for(ControlledExperiment rep : manager.getReplicates()){
			if(mgpsconfig.getDefaultBindingModel()!=null)
				bindingManager.setBindingModel(rep, mgpsconfig.getDefaultBindingModel());
			else if(rep.getExptType()!=null && rep.getExptType().getName().toLowerCase().equals("chipexo"))
				bindingManager.setBindingModel(rep, new BindingModel(BindingModel.defaultChipExoEmpiricalDistribution));
			else if(rep.getExptType()!=null && rep.getExptType().getName().toLowerCase().equals("permchipseq"))
				bindingManager.setBindingModel(rep, new BindingModelPerBase(BindingModelPerBase.defaultPermChipSeqEmpiricalDistribution));
			else
				bindingManager.setBindingModel(rep, new BindingModel(BindingModel.defaultChipSeqEmpiricalDistribution));
			repBindingModels.put(rep, new ArrayList<BindingModel>());
			repBindingModels.get(rep).add(bindingManager.getBindingModel(rep));
		}
		for(ExperimentCondition cond : manager.getConditions())
			bindingManager.updateMaxInfluenceRange(cond);
		
		//Find potential binding regions
		System.err.println("Finding potential binding regions.");
		potentialFilter = new PotentialRegionFilter(mgpsconfig, econfig, manager, bindingManager);
		List<Region> potentials = potentialFilter.execute();
		System.err.println(potentials.size()+" potential regions found.");
		if(potentials.size()==0){
			System.err.println("No potential regions - exiting.");
			System.exit(1);
		}
		potentialFilter.printPotentialRegionsToFile();
	}
	
	/**
	 * Run the mixture model to find binding events. 
	 */
	public void runMixtureModel() {
		Double[] kl;
		System.err.println("Initialzing mixture model");
		mixtureModel = new BindingMixture(gconfig, econfig, mgpsconfig, manager, bindingManager, potentialFilter);
		
		int round = 0;
		boolean converged = false;
        while (!converged){
        	
            System.err.println("\n============================ Round "+round+" ============================");
            
            //Execute the mixture model
            if(round==0)
            	mixtureModel.execute(true, true); //EM
            else
            	mixtureModel.execute(true, false); //EM
            
            //Update binding models
            String distribFilename = mgpsconfig.getOutputIntermediateDir()+File.separator+mgpsconfig.getOutBase()+"_t"+round+"_ReadDistrib";
            kl = mixtureModel.updateBindingModel(distribFilename);
            //Add new binding models to the record
            for(ControlledExperiment rep : manager.getReplicates())
    			repBindingModels.get(rep).add(bindingManager.getBindingModel(rep));
            mixtureModel.updateAlphas();
            
            //Update motifs
            mixtureModel.updateMotifs();
            
            //Update noise models
            mixtureModel.updateGlobalNoise();
            
            //Print current components
            mixtureModel.printActiveComponentsToFile();
            
            round++;
            
            //Check for convergence
            if(round>mgpsconfig.getMaxModelUpdateRounds()){
            	converged=true;
            }else{
            	converged = true;
            	for(int l=0; l<kl.length; l++)
            		converged = converged && (kl[l]<-5 || kl[l].isNaN());
            }
        }
        outFormatter.plotAllReadDistributions(repBindingModels);
        
        //ML quantification of events
        System.err.println("\n============================ ML read assignment ============================");
        mixtureModel.execute(false, false); //ML
        bindingManager.setBindingEvents(mixtureModel.getBindingEvents());
        //Update sig & noise counts in each replicate
        bindingManager.estimateSignalVsNoiseFractions(bindingManager.getBindingEvents());
        System.err.println("ML read assignment finished.");
        
        System.err.println("\n============================= Post-processing ==============================");
        
        //Statistical analysis: Enrichment over controls 
        EnrichmentSignificance tester = new EnrichmentSignificance(mgpsconfig, manager, bindingManager, mgpsconfig.getMinEventFoldChange(), econfig.getMappableGenomeLength());
		tester.execute();
        
		//Write the replicate counts to a file (needed before EdgeR differential enrichment)
		bindingManager.writeReplicateCounts(mgpsconfig.getOutputParentDir()+File.separator+mgpsconfig.getOutBase()+".replicates.counts");
		
		//Statistical analysis: inter-condition differences
		if(manager.getNumConditions()>1 && mgpsconfig.getRunDiffTests()){
			normalizer = new TMMNormalization(manager.getReplicates().size(), 0.3, 0.05);
			DifferentialEnrichment edgeR = new EdgeRDifferentialEnrichment(mgpsconfig);
			
			for(int ref=0; ref<manager.getNumConditions(); ref++){
				data = new CountsDataset(manager, bindingManager.getBindingEvents(), ref);
				//normalizer.normalize(data);
				//data.calcScMeanAndFold();
				edgeR.setFileIDname("_"+manager.getIndexedCondition(ref).getName());
				data = edgeR.execute(data);
				data.updateEvents(bindingManager.getBindingEvents(), manager);
				
				//Print MA scatters (inter-sample & inter-condition)
				//data.savePairwiseFocalSampleMAPlots(config.getOutputImagesDir()+File.separator, true);
				data.savePairwiseConditionMAPlots(mgpsconfig.getDiffPMinThres(), mgpsconfig.getOutputImagesDir()+File.separator, true);
				
				//Print XY scatters (inter-sample & inter-condition)
				data.savePairwiseFocalSampleXYPlots(mgpsconfig.getOutputImagesDir()+File.separator, true);
				data.savePairwiseConditionXYPlots(manager, bindingManager, mgpsconfig.getDiffPMinThres(), mgpsconfig.getOutputImagesDir()+File.separator, true);
			}
		}
        
        // Print final events to files
		bindingManager.writeBindingEventFiles(mgpsconfig.getOutputParentDir()+File.separator+mgpsconfig.getOutBase(), mgpsconfig.getQMinThres(), mgpsconfig.getRunDiffTests(), mgpsconfig.getDiffPMinThres());
		if(mgpsconfig.getFindingMotifs())
			bindingManager.writeMotifFile(mgpsconfig.getOutputParentDir()+File.separator+mgpsconfig.getOutBase()+".motifs");
        System.err.println("Binding event detection finished!\nBinding events are printed to files in "+mgpsconfig.getOutputParentDir()+" beginning with: "+mgpsconfig.getOutName());
        
        //Post-analysis of peaks
        EventsPostAnalysis postAnalyzer = new EventsPostAnalysis(mgpsconfig, manager, bindingManager, bindingManager.getBindingEvents(), mixtureModel.getMotifFinder());
        postAnalyzer.execute(400);
    }
	
	/**
	 * Main driver method for MultiGPS
	 * @param args
	 */
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("MultiGPS version "+MultiGPSConfig.version+"\n\n");
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		MultiGPSConfig config = new MultiGPSConfig(gcon, args);
		if(config.helpWanted()){
			System.err.println(MultiGPS.getMultiGPSArgsList());
		}else{
			
			ExperimentManager manager = new ExperimentManager(econ);
			
			//Just a test to see if we've loaded all conditions
			if(manager.getConditions().size()==0){
				System.err.println("No experiments specified. Use --expt or --design options."); System.exit(1);
			}
			
			MultiGPS gps = new MultiGPS(gcon, econ, config, manager);
			gps.runMixtureModel();
			
			manager.close();
		}
		System.err.println("Closing");
		//System.exit(0);
	}
	
	/**
	 * returns a string describing the arguments for the public version of MultiGPS. 
	 * @return String
	 */
	public static String getMultiGPSArgsList(){
		return(new String("" +
				"Copyright (C) Shaun Mahony 2012-2014\n" +
				"<http://mahonylab.org/software/multigps>\n" +
				"\n" +
				"MultiGPS comes with ABSOLUTELY NO WARRANTY.  This is free software, and you\n"+
				"are welcome to redistribute it under certain conditions.  See the MIT license \n"+
				"for details.\n"+
				"\n OPTIONS:\n" +
				" General:\n"+
				"\t--out <output file prefix>\n" +
				"\t--threads <number of threads to use>\n" +
				"\t--verbose [flag to print intermediate files and extra output]\n" +
				"\t--config <config file: all options here can be specified in a name<space>value text file, over-ridden by command-line args>\n" +
				" Genome:\n" +
				"\t--geninfo <genome info file> AND --seq <fasta seq directory reqd if using motif prior>\n" +
				" Loading Data:\n" +
				"\t--expt <file name> AND --format <SAM/BED/IDX>\n" +
				"\t--ctrl <file name (optional argument. must be same format as expt files)>\n" +
				"\t--design <experiment design file name to use instead of --expt and --ctrl; see website for format>\n"+
				"\t--fixedpb <fixed per base limit>\n" +
				"\t--poissongausspb <filter per base using a Poisson threshold parameterized by a local Gaussian sliding window>\n" +
				"\t--nonunique [flag to use non-unique reads]\n" +
				" Scaling data:\n"+
				"\t--noscaling [flag to turn off signal vs control scaling (default = scaling by regression)]\n" +
				"\t--medianscale [flag to use scaling by median (default = scaling by regression)]\n" +
				"\t--sesscale [flag to use scaling by SES (default = scaling by regression)]\n" +
				"\t--scalewin <window size for scaling procedure>\n" +
				" Running MultiGPS:\n" +
				"\t--d <binding event read distribution file>\n" +
				"\t--r <max. model update rounds>\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--minmodelupdateevents <minimum number of events to support an update>\n" +
				"\t--nomodelsmoothing [flag to turn off binding model smoothing]\n" +
				"\t--splinesmoothparam <spline smoothing parameter>\n" +
				"\t--gaussmodelsmoothing [flag to turn on Gaussian model smoothing (default = cubic spline)]\n" +
				"\t--gausssmoothparam <Gaussian smoothing std dev>\n" +
				"\t--jointinmodel [flag to allow joint events in model updates (default=do not)]\n" +
				"\t--mlconfignotshared [flag to not share component configs in the ML step]\n" +
				"\t--exclude <file of regions to ignore>\n" +
				" MultiGPS priors:\n"+
				"\t--noposprior [flag to turn off inter-experiment positional prior (default=on)]\n" +
				"\t--probshared <probability that events are shared across conditions (default=0.9)>\n" +
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--mememinw <minw arg for MEME (default=6)>\n"+
				"\t--mememaxw <maxw arg for MEME (default=18)>\n"+
				"\t--memeargs <additional args for MEME (default=  -dna -mod zoops -revcomp -nostatus)>\n"+
				" Reporting binding events:\n" +
				"\t--q <Q-value minimum (corrected p-value)>\n" +
				"\t--minfold <minimum event fold-change vs scaled control>\n" +
				"\t--nodifftests [flag to turn off differential enrichment tests]\n" +
				"\t--rpath <path to the R bin dir (default: R is in $PATH). Note that you need to install edgeR separately>\n" +
				"\t--edgerod <EdgeR overdispersion parameter>\n" +
				"\t--diffp <minimum p-value for reporting differential enrichment>\n" +
				""));
	}

}
