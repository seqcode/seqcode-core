package edu.psu.compbio.seqcode.projects.naomi;

import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;


public class ChIPexoQC {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	
	public ChIPexoQC(GenomeConfig gcon, ExptConfig econ){	
		gconfig = gcon;
		econfig = econ;
	}
	
	public void printQCMetrics(){
		
		double ncis, signalHits, noiseHits, IPstrength;
		
		ExperimentManager manager = new ExperimentManager(econfig);
				
		for(ExperimentCondition exptCond: manager.getConditions()){
			for(ControlledExperiment rep : exptCond.getReplicates()){
				ncis = rep.getSignalVsNoiseFraction();
				signalHits = rep.getSignal().getHitCount();
				noiseHits = rep.getControl().getHitCount();
				IPstrength = 1-(ncis/(signalHits/noiseHits));
				System.out.println("signal name is: "+rep.getSignal());
				System.out.println("control name is: "+rep.getControl());
				System.out.println("IP strength is: "+IPstrength);
			}
		}
	}
		
	public static void main(String[] args) {
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		ChIPexoQC exoQC = new ChIPexoQC(gconf, econf); 
		exoQC.printQCMetrics();
	}

}
