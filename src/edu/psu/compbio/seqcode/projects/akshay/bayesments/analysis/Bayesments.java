package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.BIC;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.MAPassignment;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;

/**
 * The main driver class for the Bayesments project
 * @author akshaykakumanu
 *
 */
public class Bayesments {
	
	public static void main(String[] args){
		// Build the Config class the comman line arguments
		Config c = new Config(args);
		if(c.helpWanter()){
			System.err.println("Bayesments:");
			System.err.println(c.getArgsList());
		}else{
			// Make the output directories using the provided root name
			c.makeGPSOutputDirs(c.doEMplot());
			
			ExperimentManager manager = null;
			GenomicLocations trainingData = null;
			
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
			}
			EMrunner trainer = null;
			
			double[][] bic_vals = new double[c.getMaxChromStates()-c.getMinChromStates()+1][c.getMaxFacStates()-c.getMinFacStates()+1];
			int nRows = c.getMaxChromStates()-c.getMinChromStates()+1;
			int nCols = c.getMaxFacStates()-c.getMinFacStates()+1;
			for(int i=0; i<nRows; i++){
				for(int j=0; j<nCols; j++){
					if(!c.doSimulation()){
						trainer = new EMrunner(c, trainingData, manager, i+c.getMinChromStates(), j+c.getMinFacStates(), false);
					}else{
						trainer = new EMrunner(c, i+c.getMinChromStates(), j+c.getMinFacStates(), false);
					}
					trainer.trainModel();
					BIC temp_bic = new BIC(c, trainer.getModel(), i+c.getMinChromStates(), j+c.getMinFacStates());
					bic_vals[i][j] = temp_bic.calculateBicScore();
				}
			}
			
			int nChromStates = BayesmentsSandbox.getMinIndex(bic_vals).car()+c.getMinChromStates();
			int nFacStates = BayesmentsSandbox.getMinIndex(bic_vals).cdr()+c.getMinFacStates();
			
			
			trainer = null;
			
			
			if(!c.doSimulation()){
				trainer = new EMrunner(c, trainingData, manager, nChromStates, nFacStates,c.doRegularization());
				
			}else{
				trainer = new EMrunner(c, nChromStates, nFacStates,c.doRegularization());
			}
			//Initialize EMrunner and train the model
			trainer.trainModel();
			
			//Do MAP assignament
			if(!c.doSimulation()){
				MAPassignment map =  new MAPassignment(trainer.getModel(), c, trainingData.getLocations(), c.getNumChrmStates(), c.getNumFacStates());
				map.execute(true);
			}else{
				MAPassignment map =  new MAPassignment(trainer.getModel(), c, null, c.getNumChrmStates(), c.getNumFacStates());
				map.execute(false);
			}
			
			
			
			// Print all the learned parameters
			//System.out.println("PI-C values\n");
			//BayesmentsSandbox.printArray(trainer.getPIj(), "chrom_state");
		
			System.out.println("MU-C values\n");
			BayesmentsSandbox.printArray(trainer.getMUc(),"MUc" , "MUc", trainer.getModel().getConditionNames());
		
			System.out.println("MU-F values\n");
			BayesmentsSandbox.printArray(trainer.getMUf(),"MUf" , "MUf", trainer.getModel().getConditionNames());
			if(!c.runOnlyChrom()){
				System.out.println("MU-S values\n");
				BayesmentsSandbox.printArray(trainer.getMUs(),"MUS" , "MUS", trainer.getModel().getConditionNames());
			}
		
			System.out.println("SIGMA-C values\n");
			BayesmentsSandbox.printArray(trainer.getSIGMAc(),"SIGMAc" , "SIGMAc", trainer.getModel().getConditionNames());
		
			System.out.println("SIGMA-F values\n");
			BayesmentsSandbox.printArray(trainer.getSIGMAf(),"SIGMAf" , "SIGMAf", trainer.getModel().getConditionNames());
			
			if(!c.runOnlyChrom()){
				System.out.println("SIGMA-S values\n");
				BayesmentsSandbox.printArray(trainer.getSIGMAs(),"SIGMAs" , "SIGMAs", trainer.getModel().getConditionNames());
			}
		
			System.out.println("Bjk values\n");
			BayesmentsSandbox.printArray(trainer.getBjk(),"chromatin_State" , "factor_state", trainer.getModel().getConditionNames());
			
			if(c.doRegularization()){
				System.out.println("Chromatin Weight values\n");
				BayesmentsSandbox.printArray(trainer.getChromWeights(), "Chromatin Feature");
				
				if(!c.runOnlyChrom()){
					System.out.println("Sequence Weight values\n");
					BayesmentsSandbox.printArray(trainer.getSeqWeights(), "Sequence Feature");
				}
				
			}
		}
	}	
}
