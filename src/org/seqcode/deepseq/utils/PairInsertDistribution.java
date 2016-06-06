package org.seqcode.deepseq.utils;

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;
import org.seqcode.utils.NotFoundException;
import org.seqcode.utils.RealValuedHistogram;


/**
 * Utility to print a histogram of paired read inferred fragment sizes
 * 
 *  Unit testing:
 *  	
 *  
 * @author mahony
 *
 */
public class PairInsertDistribution {

	private GenomeConfig gconfig;
	private ExptConfig econfig;
	private ExperimentManager manager=null;
	private HashMap<Sample, RealValuedHistogram> histograms;
	private int histoMax=2000;
	private int histoBins = 200;
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if(args.length==0 || gcon.helpWanted()){
			System.err.println("PairInsertDistribution\n"+gcon.getArgsList()+"\n"+econ.getArgsList());
			System.exit(1);
		}
		
        PairInsertDistribution pid = new PairInsertDistribution(gcon, econ);
        pid.execute();
        
        pid.close();
	}
	
	public PairInsertDistribution(GenomeConfig g, ExptConfig e){
		gconfig = g;
		econfig = e;
		econfig.setLoadPairs(true);
		histograms = new HashMap<Sample, RealValuedHistogram>();
		
		manager = new ExperimentManager(this.econfig);
		for(Sample samp : manager.getSamples())
			histograms.put(samp, new RealValuedHistogram(0, histoMax, histoBins));
	}
	
	public void execute(){
		Genome gen = gconfig.getGenome();

		Iterator<Region> testRegions = new ChromosomeGenerator().execute(gen);
		while (testRegions.hasNext()) {
			Region currReg = testRegions.next();
			if(!currReg.getChrom().endsWith("_random")){
				System.err.println(currReg.getLocationString());
				for(Sample samp : manager.getSamples()){
					List<StrandedPair> sps = samp.getPairs(currReg);
					for(StrandedPair pair : sps){
						int fs = pair.getFragmentSize();
						if(fs !=-1){
							histograms.get(samp).addValue(fs);
						}
					}
				}
			}
		}
		for(Sample samp : manager.getSamples()){
			System.out.println(samp.getName()+"\nTotalPairs = "+samp.getPairCount()+", UniquePairs = "+samp.getUniquePairCount()+"\nBin\tCount");
			histograms.get(samp).printContents();
		}
	}
	
	public void close(){
		if(manager!=null)
			manager.close();
	}
}
