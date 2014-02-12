package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

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
			
			// Store the reads by building the ExperimentManager
			System.out.println("Loading Reads\n");
			ExperimentManager manager = new ExperimentManager(c,c.loadReads());
			
			// Read and store the training data in a GenomicLocations object
			System.out.println("Filling Training Data\n");
			GenomicLocations trainingData = new GenomicLocations(manager,c);
			
			// Plot the cumulative plots of the training data
			System.out.println("Plotting the traning data");
			trainingData.plotData(c, manager);
			
			// Initialize the EM training object
			System.out.println("Initializing EM\n");
			EMtrain EM = new EMtrain(c,trainingData,manager);
			
			//Run EM
			System.out.println("Running EM\n");
			//EM.runEM();
			
			// Perform MAP assignment
			MAPassignment assignment = new MAPassignment(EM, c);
			assignment.execute();
			
			// Print all the learned parameters
			System.out.println("PI-C values\n");
			BayesmentsSandbox.printArray(EM.getPIj(), "chrom_state");
		
			System.out.println("MU-C values\n");
			BayesmentsSandbox.printArray(EM.getMUc(),"MUc" , "MUc", manager);
		
			System.out.println("MU-F values\n");
			BayesmentsSandbox.printArray(EM.getMUf(),"MUf" , "MUf", manager);
		
			System.out.println("SIGMA-C values\n");
			BayesmentsSandbox.printArray(EM.getSIGMAc(),"SIGMAc" , "SIGMAc", manager);
		
			System.out.println("SIGMA-F values\n");
			BayesmentsSandbox.printArray(EM.getSIGMAf(),"SIGMAf" , "SIGMAf", manager);
		
			System.out.println("Bjk values\n");
			BayesmentsSandbox.printArray(EM.getBjk(),"chromatin_State" , "factor_state", manager);
		}
	}	
}
