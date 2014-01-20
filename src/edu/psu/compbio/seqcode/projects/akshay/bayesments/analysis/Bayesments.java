package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;

public class Bayesments {
	
	public static void main(String[] args){
		
		Config c = new Config(args);
		if(c.helpWanter()){
			System.err.println("MultiGPS:");
			System.err.println(c.getArgsList());
		}else{
			c.makeGPSOutputDirs(c.doEMplot());
			System.out.println("Loading Reads\n");
			ExperimentManager manager = new ExperimentManager(c,c.loadReads());
			System.out.println("Filling Training Dara\n");
			GenomicLocations trainingData = new GenomicLocations(manager,c);
			
			System.out.println("Plotting the traning data");
			trainingData.plotData(c, manager);
			
			System.out.println("Initializing EM\n");
			EMtrain EM = new EMtrain(c,trainingData);
			
			System.out.println("Running EM\n");
			EM.runEM();
		
			System.out.println("PI-C values\n");
			BayesmentsSandbox.printArray(EM.getPIj(), "chrom_state");
		
			System.out.println("MU-C values\n");
			BayesmentsSandbox.printArray(EM.getMUc(),"chrom_State" , "condition");
		
			System.out.println("MU-F values\n");
			BayesmentsSandbox.printArray(EM.getMUf(),"factor_State" , "condition");
		
			System.out.println("SIGMA-C values\n");
			BayesmentsSandbox.printArray(EM.getSIGMAc(),"chrom_State" , "condition");
		
			System.out.println("SIGMA-F values\n");
			BayesmentsSandbox.printArray(EM.getSIGMAf(),"factor_State" , "condition");
		
			System.out.println("Bjk values\n");
			BayesmentsSandbox.printArray(EM.getBjk(),"chromatin_State" , "factor_state");
		}
	}	
}
