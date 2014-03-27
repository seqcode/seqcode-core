package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.BIC;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.MAPassignment;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.Sequences;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;
import edu.psu.compbio.seqcode.projects.shaun.EventMetaMaker;

/**
 * The main driver class for the Bayesments project
 * @author akshaykakumanu
 *
 */
public class Bayesments {
	
	public static void main(String[] args) throws IOException{
		// Build the Config class the comman line arguments
		Config c = new Config(args);
		if(c.helpWanter()){
			System.err.println("Bayesments:");
			System.err.println(c.getArgsList());
		}else{
			// Make the output directories using the provided root name
			c.makeGPSOutputDirs(true);
			
			ExperimentManager manager = null;
			GenomicLocations trainingData = null;
			Sequences seqs = null;
			List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
			
			//Loading/reading tags
			if(!c.doSimulation()){
				// Store the reads by building the ExperimentManager
				System.out.println("Loading Reads\n");
				manager = new ExperimentManager(c,c.loadReads());
				
				// Read and store the training data in a GenomicLocations object
				System.out.println("Filling Training Data\n");
				trainingData = new GenomicLocations(manager,c);
				
				// Plot the cumulative plots of the training data
				System.out.println("Plotting the traning data");
				trainingData.plotData(c, manager);
				
				if(!c.runOnlyChrom()){
					//Reading the genome sequence and intializing the sequences class object
					System.out.println("Reading the genome sequence");
					seqs = new Sequences(c,trainingData.getLocations());
					List<File> top_locs = c.getTopLocations();
					for(int f=0; f<top_locs.size(); f++){
						List<Point> top_regions = EventMetaMaker.loadPoints(top_locs.get(f), c.getGenome());
						MotifPlatform motif_profiler = new MotifPlatform(c,top_regions);
						motif_profiler.findMotifs();
						List<WeightMatrix> temp_motifs = motif_profiler.getMotifs();
						for(int m=0; m< temp_motifs.size(); m++){
							motifs.add(temp_motifs.get(m));
						}
					}
					seqs.setMotifs(motifs);
					seqs.setXc();
					seqs.plotSeqScores();
				}
			}
			
			
			EMtrain trainer = null;
			
			double[][] bic_vals = new double[c.getMaxChromStates()-c.getMinChromStates()+1][c.getMaxFacStates()-c.getMinFacStates()+1];
			double[][] aic_vals = new double[c.getMaxChromStates()-c.getMinChromStates()+1][c.getMaxFacStates()-c.getMinFacStates()+1];
			int nRows = c.getMaxChromStates()-c.getMinChromStates()+1;
			int nCols = c.getMaxFacStates()-c.getMinFacStates()+1;
			for(int i=0; i<nRows; i++){
				for(int j=0; j<nCols; j++){
					if(!c.doSimulation()){
						if(c.runOnlyChrom() || motifs == null){
							trainer = new EMtrain(c, trainingData, null, manager,i+c.getMinChromStates(), j+c.getMinFacStates(), c.doRegularization(), false, false);
						}else{
							trainer = new EMtrain(c, trainingData, seqs, manager, i+c.getMinChromStates(), j+c.getMinFacStates(), c.doRegularization(), false, true);
						}
					}else{
						trainer = new EMtrain(c, i+c.getMinChromStates(), j+c.getMinFacStates(), c.doRegularization(), false);
					}
					trainer.runEM(false);
					BIC temp_bic = new BIC(c, trainer, i+c.getMinChromStates(), j+c.getMinFacStates());
					bic_vals[i][j] = temp_bic.calculateBicScore();
					aic_vals[i][j] = temp_bic.calculateAicScore();
				}
			}
			
			int nChromStates = BayesmentsSandbox.getMinIndex(bic_vals).car()+c.getMinChromStates();
			int nFacStates = BayesmentsSandbox.getMinIndex(bic_vals).cdr()+c.getMinFacStates();
			
			//Printing 
			
			System.out.println("-------------------------------------------------------------BIC values:----------------------------------------------------------");
			BayesmentsSandbox.printArray(bic_vals,"chromatin_State" , "factor_state", trainer.getConditionNames());
			System.out.println("----------------------------------------------------------------------------------------------------------------------------------");
			
			System.out.println("-------------------------------------------------------------AIC values:----------------------------------------------------------");
			BayesmentsSandbox.printArray(aic_vals,"chromatin_State" , "factor_state", trainer.getConditionNames());
			System.out.println("----------------------------------------------------------------------------------------------------------------------------------");
			
			trainer = null;
			
			
			if(!c.doSimulation()){
				if(c.runOnlyChrom() || motifs == null){
					trainer = new EMtrain(c, trainingData, null, manager,nChromStates,nFacStates,c.doRegularization(), true, false);
				}else{
					trainer = new EMtrain(c, trainingData, seqs, manager,nChromStates,nFacStates, c.doRegularization(), true, true);
				}
				
			}else{
				trainer = new EMtrain(c, nChromStates, nFacStates, c.doRegularization(), true);
			}
			//Initialize EMrunner and train the model
			trainer.runEM(c.doEMplot());
			
			//Do MAP assignament
			if(!c.doSimulation()){
				MAPassignment map =  new MAPassignment(trainer, c, trainingData.getLocations(), nChromStates, nFacStates);
				map.execute(true);
			}else{
				MAPassignment map =  new MAPassignment(trainer, c, null, nChromStates, nFacStates);
				map.execute(false);
			}
			
			System.out.println("MU-C values\n");
			BayesmentsSandbox.printArray(trainer.getMUc(),"MUc" , "MUc", trainer.getConditionNames());
		
			System.out.println("MU-F values\n");
			BayesmentsSandbox.printArray(trainer.getMUf(),"MUf" , "MUf", trainer.getConditionNames());
			if(!c.runOnlyChrom() && motifs != null){
				System.out.println("MU-S values\n");
				BayesmentsSandbox.printArray(trainer.getMUs(),"MUS" , "MUS", trainer.getConditionNames());
			}
		
		
			System.out.println("Bjk values\n");
			BayesmentsSandbox.printArray(trainer.getBjk(),"chromatin_State" , "factor_state", trainer.getConditionNames());
			
			if(c.doRegularization()){
				System.out.println("Chromatin Weight values\n");
				BayesmentsSandbox.printArray(trainer.getChromWeights(), "Chromatin Feature");
				
				if(!c.runOnlyChrom() && motifs != null){
					System.out.println("Sequence Weight values\n");
					BayesmentsSandbox.printArray(trainer.getSeqWeights(), "Sequence Feature");
				}
				
			}
			// Merging chromatin and sequence matrices
			if(trainer.getSeqStateStatus()){
				double[][] chromatin_Mu = new double[trainer.getMUc().length][trainer.getMUc()[0].length+trainer.getMUs()[0].length];
				for(int row=0; row<chromatin_Mu.length;row++){
					for(int col=0; col<trainer.getMUc()[0].length;col++){
						chromatin_Mu[row][col] = trainer.getMUc()[row][col];
					}
				}
				for(int row = 0; row< chromatin_Mu.length; row++){
					for(int col=trainer.getMUc()[0].length;col<chromatin_Mu[0].length; col++){
						chromatin_Mu[row][col] = trainer.getMUs()[row][col-trainer.getMUc()[0].length];
					}
				}
			
				// Generating xlabs 
			
				String[] chromatin_name = new String[chromatin_Mu[0].length];
				for(int col=0; col<trainer.getMUc()[0].length; col++ ){
					chromatin_name[col] = trainer.getConditionNames()[col];
				}
				for(int col= trainer.getMUc()[0].length; col<chromatin_Mu[0].length; col++){
					chromatin_name[col] = "Motif "+Integer.toString(col-trainer.getMUc()[0].length);
				}
			
				//plotting MUc and MUs heats
			
				HeatMapper map = new HeatMapper(c, chromatin_Mu, "Experimental_Track", "State", chromatin_name, "Mu_chromatin");
				map.plot(new Color(221,20,20),true);
			
				// Plotting factor states
			}else{
				HeatMapper map = new HeatMapper(c, trainer.getMUc(), "Experimental_Track", "State", trainer.getConditionNames(), "Mu_chromatin");
			}
			
			HeatMapper map = new HeatMapper(c, trainer.getMUf(), "Factor-Expt", "State", null, "Mu_factor");
			map.plot(new Color(20,221,20),true);
			
			//Plotting transitions
			
			map =  new HeatMapper(c, trainer.getBjk(), "Fac-State", "Chrom-State", null, "Transitions");
			map.plot(new Color(20,20,221),false);
			
		}
	}	
}
