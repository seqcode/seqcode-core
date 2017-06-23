package org.seqcode.tutorials;

import java.util.Iterator;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;


/**
 * Simple package testing class to cache a set of experiments and print total hit counts.
 * 
 *  Unit testing:
 *  	java -Xmx2G org.seqcode.tutorials.SeqExperimentLoadingExample --species "Mus musculus;mm9" --rdbexptC1 "ES2MN Day0(iCdx2.V5+Dox) iCdx2.V5 Ainv15_iCdx2.V5;1;bowtie_unique" --fixedpb 100000 (--nocache)
 *  	Result: C1:rep1	4700200.0	4700200.0
 *  
 * @author mahony
 *
 */
public class SeqExperimentLoadingExample {

	final int MAXSECTION = 50000000;
	GenomeConfig gconfig;
	ExptConfig econfig;
	ExperimentManager manager;
	
	public static void main(String[] args){
		//GenomeConfig and ExptConfig read various command-line options in order to load requested genome and sequencing datasets (resp.)
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		if(args.length==0 || gconfig.helpWanted()){
			System.err.println("TestExptLoading\n"+gconfig.getArgsList()+"\n"+econfig.getArgsList());
		}
		
		//Everything required to get the data requested from the command-line args is now in the *Config objects
		SeqExperimentLoadingExample tel = new SeqExperimentLoadingExample(gconfig, econfig);
		tel.execute();
		
		//Always close connections at the end of processing
		tel.close();
	}
	public SeqExperimentLoadingExample(GenomeConfig g, ExptConfig e){
		this.gconfig = g;
		this.econfig = e;
	}
	
	/**
	 * The functionality here is trivial - just counts the number of hits in the experiment.
	 * The point is to demonstrate iteration over chromosomes and
	 * to test read loading per chromosome. 
	 */
	public void execute(){
		//ExperimentManager initializes the data structures associated with the experimental data pointed to on cmd line
		manager = new ExperimentManager(econfig);
		//Appropriate Genome will already be loaded by GenomeConfig
		Genome gen = gconfig.getGenome();
		
		
		//If we have multiple experiments, process one at a time
		for(ControlledExperiment expt : manager.getReplicates()){
			float totalCount = 0;
			//Iterate over chromosomes in the genome
			Iterator<Region> chroms = new ChromosomeGenerator().execute(gen);
			while (chroms.hasNext()) {
				double chrCount =0;
				Region currentRegion = chroms.next();
				//Split the chromosome up into large chunks
	            for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION){
	                int y = x+MAXSECTION; 
	                if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
	                Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
	                //Get hits in this subregion from the current experiment's signal track
					List<StrandedBaseCount> ipHits = expt.getSignal().getBases(currSubRegion);
					for(StrandedBaseCount z : ipHits){
						totalCount+=z.getCount();
						chrCount += z.getCount();
					}
	            }
	            System.out.println("\t"+currentRegion.getLocationString()+"\t"+chrCount);
			}
			System.out.println(expt.getName()+"\t"+expt.getSignal().getHitCount()+"\t"+totalCount);
			
		}
	}
	
	protected void close(){
		//Always close database connections (if applicable) via a call to ExperimentManager.close()
		manager.close();
	}

}
