package org.seqcode.deepseq.utils;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;


/**
 * Simple package testing class to cache a set of experiments and print total hit counts.
 * 
 *  Unit testing:
 *  	java -Xmx2G org.seqcode.deepseq.utils.TestExptLoading --species "Mus musculus;mm9" --rdbexptC1 "ES2MN Day0(iCdx2.V5+Dox) iCdx2.V5 Ainv15_iCdx2.V5;1;bowtie_unique" --fixedpb 100000 (--nocache)
 *  	Result: C1:rep1	4700200.0	4700200.0
 *  
 * @author mahony
 *
 */
public class TestExptLoading {

	final int MAXSECTION = 50000000;
	GenomeConfig gconfig;
	ExptConfig econfig;
	ExperimentManager manager;
	
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		if(args.length==0 || gconfig.helpWanted()){
			System.err.println("TestExptLoading\n"+gconfig.getArgsList()+"\n"+econfig.getArgsList());
		}
		
		TestExptLoading tel = new TestExptLoading(gconfig, econfig);
		tel.execute();
		
		tel.close();
	}
	public TestExptLoading(GenomeConfig g, ExptConfig e){
		this.gconfig = g;
		this.econfig = e;
	}
	
	/**
	 * The functionality here is trivial - the point is to demonstrate iteration over chromosomes and
	 * to test read loading per chromosome. 
	 */
	public void execute(){
		manager = new ExperimentManager(econfig);
		Genome gen = gconfig.getGenome();
		
		Map<ControlledExperiment, Double> meanFragmentCoverage = new HashMap<ControlledExperiment, Double>();
		//If we have multiple experiments, process one at a time
		for(ControlledExperiment expt : manager.getReplicates()){
			float totalCount = 0;
			Iterator<Region> chroms = new ChromosomeGenerator().execute(gen);
			while (chroms.hasNext()) {
				double chrCount =0;
				Region currentRegion = chroms.next();
				//Split the job up into large chunks
	            for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION){
	                int y = x+MAXSECTION; 
	                if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
	                Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
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
		manager.close();
	}
	
	protected void close(){
		manager.close();
	}

}
