package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class Bayesments {
	
	public static void run(String[] args){
		Config c = new Config(args);
		c.makeGPSOutputDirs(c.doEMplot());
		
		ExperimentManager manager = new ExperimentManager(c,c.loadReads());
		GenomicLocations trainingData = new GenomicLocations(manager,c);
		
		EMtrain EM = new EMtrain(c,trainingData);
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
