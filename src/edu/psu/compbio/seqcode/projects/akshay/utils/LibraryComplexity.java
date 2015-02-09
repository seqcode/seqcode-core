package edu.psu.compbio.seqcode.projects.akshay.utils;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;

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
	
		LibraryComplexity runner = new LibraryComplexity(econ,gcon);
		runner.execute();
		
	}
	/**
	 * 
	 */
	public void execute(){
		ExperimentManager manager = new ExperimentManager(econfig);
		try{
			if(manager.getSamples().size()>0){
				Sample sample = manager.getSamples().get(0);
				double totHits = sample.getHitCount();
				double uniqPos = sample.getHitPositionCount();
				double complexity = uniqPos/totHits;
				System.out.println("Total hits: "+Double.toString(totHits));
				System.out.println("Total Unique positions: "+Double.toString(uniqPos));
				System.out.println("Complexity: "+ Double.toString(complexity));
			}else{
				throw new Exception("Provide a Readdb experiment");
				
			}
			
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	

}
