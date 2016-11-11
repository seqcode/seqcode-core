package org.seqcode.projects.galaxyexo;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.ml.math.stats.StatUtil;

/**
 * A method to quantify reads enrichment around plus one nucleosome relative to control experiment. 
 * 
 * input : reference point
 *       : signal and control experiment
 * 
 * @author naomi yamada
 */

public class NucleosomeEnrichmentProfiler {
	protected ExperimentManager manager;
	protected ExptConfig exptConfig;
	protected FeatureCountsLoader featureCountsLoader;
	protected String outbase;
	protected int smoothingWidth;
	
	public NucleosomeEnrichmentProfiler(String base,ExptConfig econf, ExperimentManager man, FeatureCountsLoader fcloader, int smoothing){
		outbase = base;
		featureCountsLoader = fcloader;
		exptConfig = econf;
		manager = man;		
		smoothingWidth=smoothing;
	}
	
	public void execute() throws FileNotFoundException{
		
		File outFile = new File(outbase);
		PrintWriter writer = new PrintWriter(outFile);
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){				
				if ( !rep.hasControl()){
					System.err.println("Please provide the control experiment to get statistics.");
					System.exit(0);
				}				
				double[] composite = featureCountsLoader.sampleComposite(rep); 
				double[] contComposite = featureCountsLoader.controlComposite();
				
				// smooth composite plot
				double[] smoothedSignal = StatUtil.gaussianSmoother(composite, smoothingWidth);
				double[] smoothedControl = StatUtil.gaussianSmoother(contComposite, smoothingWidth);
				
				double ratio = getMaxCounts(smoothedSignal)*getMinCounts(smoothedControl)/(getMinCounts(smoothedSignal)*getMaxCounts(smoothedControl));
				
				writer.println("sample "+rep.getSignal().getName()+"\tratio: "+ratio);
			}
		}
		writer.close();
	}
	
	public double getMaxCounts(double[] array){
		double maxCounts = 0;
		for (int i = 0; i < array.length; i++){
			if (array[i] > maxCounts)
				maxCounts = array[i];
		}
		return maxCounts;
	}
	
	public double getMinCounts(double[] array){
		double minCounts = Double.MAX_VALUE;
		for (int i = 0; i < array.length; i++){
			if (array[i] < minCounts)
				minCounts = array[i];
		}
		return minCounts;
	}	
	
	public static void main(String[] args) throws FileNotFoundException{
		ArgParser ap = new ArgParser(args);		
        if((!ap.hasKey("peaks") && !ap.hasKey("regions")) ) { 
        	System.err.println("please input peak files and region files.");
            System.err.println("Usage:\n " +
                               "NucleosomeEnrichmentProfiler\n " +
                               "--species <organism;genome> OR\n" +
                               "--geninfo <genome info> AND --seq <path to seqs>\n" +
                               "--peaks <file containing coordinates of peaks> \n" +
                               "\nOPTIONS:\n" +
                               "--out <output file path+name or name> \n" +
                               "--win <window of sequence to take around peaks> \n" +
                               "--bai <file path to bai file> \n" +
                               "--readshift \n" +
                               "--smooth \n" +
                               "");
            System.exit(0);
        }
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		
		// parsing command line arguments	
		int win = Args.parseInteger(args, "win", 500);
		int fivePrimeShift = Args.parseInteger(args,"readshift", 7);
		int smooth = Args.parseInteger(args,"smooth", 10);
		List<StrandedPoint> spoints = RegionFileUtilities.loadStrandedPointsFromFile(gconf.getGenome(), ap.getKeyValue("peaks"));
		// Get outdir and outbase and make them;
		String outbase = Args.parseString(args, "out", "nuc_enrichment.txt");
		
		FeatureCountsLoader fcLoader = new FeatureCountsLoader(gconf, spoints, win);
		fcLoader.setFivePrimeShift(fivePrimeShift);

		NucleosomeEnrichmentProfiler profile = new NucleosomeEnrichmentProfiler(outbase,econf,manager, fcLoader,smooth); 	
		profile.execute();
		
		manager.close();
	}
}
