package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class Bayesments {
	
	public static void run(String[] args){
		
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
			
			System.out.println("Initializing EM\n");
			EMtrain EM = new EMtrain(c,trainingData);
			
			System.out.println("Running EM\n");
			EM.runEM();
		
			System.out.println("PI-C values\n");
			System.out.println(EM.getPIj());
		
			System.out.println("MU-C values\n");
			System.out.println(EM.getMUc());
		
			System.out.println("MU-F values\n");
			System.out.println(EM.getMUf());
		
			System.out.println("SIGMA-C values\n");
			System.out.println(EM.getSIGMAc());
		
			System.out.println("SIGMA-F values\n");
			System.out.println(EM.getSIGMAf());
		
			System.out.println("Bjk values\n");
			System.out.println(EM.getBjk());
		}
	}	
}
