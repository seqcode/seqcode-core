package edu.psu.compbio.seqcode.projects.multigps.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.EnrichmentSignificance;
import edu.psu.compbio.seqcode.projects.multigps.framework.OutputFormatter;
import edu.psu.compbio.seqcode.projects.multigps.framework.PotentialRegionFilter;
import edu.psu.compbio.seqcode.projects.multigps.mixturemodel.BindingMixture;
import edu.psu.compbio.seqcode.projects.multigps.stats.CountsDataset;
import edu.psu.compbio.seqcode.projects.multigps.stats.DESeqDifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.DifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.EdgeRDifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.Normalization;
import edu.psu.compbio.seqcode.projects.multigps.stats.TMMNormalization;
import edu.psu.compbio.seqcode.projects.multigps.utilities.EventsPostAnalysis;

public class MultiGPS {

	protected ExperimentManager manager;
	protected Config config;
	protected BindingMixture mixtureModel;
	protected PotentialRegionFilter potentialFilter;
	protected OutputFormatter outFormatter;
	private CountsDataset data;
	protected Normalization normalizer;
	protected Map<ControlledExperiment, List<BindingModel>> repBindingModels;
	
	public MultiGPS(Config c, ExperimentManager eMan){
		manager = eMan;
		config = c;
		config.makeGPSOutputDirs(true);
		outFormatter = new OutputFormatter(config);
		
		//Initialize the binding model record
		repBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();
		for(ControlledExperiment rep : manager.getExperimentSet().getReplicates()){
			repBindingModels.put(rep, new ArrayList<BindingModel>());
			repBindingModels.get(rep).add(rep.getBindingModel());
		}
		
		//Find potential binding regions
		System.err.println("Finding potential binding regions.");
		potentialFilter = new PotentialRegionFilter(config, manager);
		List<Region> potentials = potentialFilter.execute();
		System.err.println(potentials.size()+" potential regions found.");
		potentialFilter.printPotentialRegionsToFile();
	}
	
	/**
	 * Run the mixture model to find binding events. 
	 */
	public void runMixtureModel() {
		Double[] kl;
		System.err.println("Initialzing mixture model");
		mixtureModel = new BindingMixture(config, manager, potentialFilter);
		
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
            String distribFilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_t"+round+"_ReadDistrib";
            kl = mixtureModel.updateBindingModel(distribFilename);
            //Add new binding models to the record
            for(ControlledExperiment rep : manager.getExperimentSet().getReplicates())
    			repBindingModels.get(rep).add(rep.getBindingModel());
            mixtureModel.updateAlphas();
            
            //Update motifs
            mixtureModel.updateMotifs();
            
            //Update noise models
            mixtureModel.updateGlobalNoise();
            
            //Print current components
            mixtureModel.printActiveComponentsToFile();
            
            round++;
            
            //Check for convergence
            if(round>config.getMaxModelUpdateRounds()){
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
        manager.setEvents(mixtureModel.getBindingEvents());
        //Update sig & noise counts in each replicate
        manager.estimateSignalProportion(manager.getEvents());
        System.err.println("ML read assignment finished.");
        
        System.err.println("\n============================= Post-processing ==============================");
        
        //Statistical analysis: Enrichment over controls 
        EnrichmentSignificance tester = new EnrichmentSignificance(config, manager.getExperimentSet(), manager.getEvents(), config.getMinEventFoldChange(), config.getMappableGenomeLength());
		tester.execute();
        
		//Write the replicate counts to a file (needed before EdgeR differential enrichment)
		manager.writeReplicateCounts(config.getOutputParentDir()+File.separator+config.getOutBase()+".replicates.counts");
		
		//Statistical analysis: inter-condition differences
		if(manager.getNumConditions()>1 && config.getRunDiffTests()){
			normalizer = new TMMNormalization(manager.getExperimentSet().getReplicates().size(), 0.3, 0.05);
			DifferentialEnrichment edgeR = new EdgeRDifferentialEnrichment(config);
			
			for(int ref=0; ref<manager.getNumConditions(); ref++){
				data = new CountsDataset(manager, manager.getEvents(), ref);
				//normalizer.normalize(data);
				//data.calcScMeanAndFold();
				data = edgeR.execute(data);
				data.updateEvents(manager.getEvents(), manager);
				
				//Print MA scatters (inter-sample & inter-condition)
				//data.savePairwiseFocalSampleMAPlots(config.getOutputImagesDir()+File.separator, true);
				data.savePairwiseConditionMAPlots(config.getDiffPMinThres(), config.getOutputImagesDir()+File.separator, true);
				
				//Print XY scatters (inter-sample & inter-condition)
				data.savePairwiseFocalSampleXYPlots(config.getOutputImagesDir()+File.separator, true);
				data.savePairwiseConditionXYPlots(manager, config.getDiffPMinThres(), config.getOutputImagesDir()+File.separator, true);
			}
		}
        
        // Print final events to files
		manager.writeBindingEventFiles(config.getOutputParentDir()+File.separator+config.getOutBase());
		manager.writeMotifFile(config.getOutputParentDir()+File.separator+config.getOutBase()+".motifs");
        System.err.println("Finished! Binding events are printed to files in "+config.getOutputParentDir()+" beginning with: "+config.getOutName());
        
        //Post-analysis of peaks
        EventsPostAnalysis postAnalyzer = new EventsPostAnalysis(config, manager, manager.getEvents(), mixtureModel.getMotifFinder());
        postAnalyzer.execute(400);
    }
	
	/**
	 * Main driver method for MultiGPS
	 * @param args
	 */
	public static void main(String[] args){
		
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("MultiGPS:");
			System.err.println(config.getArgsList());
		}else{
			ExperimentManager manager = new ExperimentManager(config);
			
			//System.err.println("Welcome to MultiGPS:\n\tStarted at "+new Date());
			//Just a test to see if we've loaded all conditions
			ExperimentSet eset = manager.getExperimentSet();
			if(eset.getConditions().size()==0){
				System.err.println("No experiments specified."); System.exit(1);
			}
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
			
			MultiGPS gps = new MultiGPS(config, manager);
			gps.runMixtureModel();
			
			manager.close();
		}
	}
}
