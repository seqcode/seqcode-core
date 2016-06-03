package org.seqcode.projects.akshay.utils;

import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.GenomeConfig;

public class LibraryComplexity {
	
	
	private ExptConfig econfig;
	private GenomeConfig gconfig;
	
	
	public LibraryComplexity(ExptConfig econ, GenomeConfig gcon) {
		econfig =econ;
		gconfig = gcon;
	}
	
	
	public static void main(String args[]){
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(),args);
		econ.setPerBaseReadFiltering(false);
	
		LibraryComplexity runner = new LibraryComplexity(econ,gcon);
		runner.execute();
		
		
	}
	/**
	 * 
	 */
	public void execute(){
		ExperimentManager manager = new ExperimentManager(econfig);
		
		if(manager.getSamples().size()>0){
			Sample sample = manager.getSamples().get(0);
			double totHits = sample.getHitCount();
			double uniqPos = sample.getHitPositionCount();
			double complexity = uniqPos/totHits;
			System.out.printf("Total hits: %.0f\n",totHits);
			System.out.printf("Total Unique positions: %.0f\n",uniqPos);
			System.out.printf("Complexity: %f\n",complexity);
			manager.close();
		}else{
			try{
				throw new Exception("Provide a Readdb experiment");
			}catch(Exception e){
				e.printStackTrace();
			}
		}
			
			
			
	}


}
