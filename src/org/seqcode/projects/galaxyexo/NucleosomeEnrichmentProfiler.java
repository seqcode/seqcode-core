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
import org.seqcode.math.stats.StatUtil;

/**
 * Utility to quantify a nucleosome "valley" of a signal experiment relative to
 * a control experiment. It creates a composite plot using reference points
 * (such as plus one nucleosome positions), and finds the maximum and minimum
 * read heights. The fold difference in the maximum and minimum read heights are
 * compared to the control experiment to get the relative fold differences in
 * read height.
 * 
 * Input: - Genome - Signal experiment - Control experiment - Peak locations
 * (plus one nucleosome) - Window around peaks in which to count reads Output: -
 * A text file containing the enrichmemt
 * 
 * @author naomi yamada
 */

public class NucleosomeEnrichmentProfiler {
	protected ExperimentManager manager;
	protected ExptConfig exptConfig;
	protected FeatureCountsLoader featureCountsLoader;
	protected String outbase;
	protected int smoothingWidth;

	public NucleosomeEnrichmentProfiler(String base, ExptConfig econf, ExperimentManager man,
			FeatureCountsLoader fcloader, int smoothing) {
		outbase = base;
		featureCountsLoader = fcloader;
		exptConfig = econf;
		manager = man;
		smoothingWidth = smoothing;
	}

	public void execute() throws FileNotFoundException {

		File outFile = new File(outbase + File.separator + "nucleosome_enrichment.txt");
		outFile.getParentFile().mkdirs();
		PrintWriter writer = new PrintWriter(outFile);

		for (ExperimentCondition condition : manager.getConditions()) {
			for (ControlledExperiment rep : condition.getReplicates()) {
				if (!rep.hasControl()) {
					System.err.println("Please provide the control experiment to get statistics.");
					System.exit(0);
				}
				double[] composite = featureCountsLoader.sampleComposite(rep);
				double[] contComposite = featureCountsLoader.controlComposite();

				// smooth composite plot
				double[] smoothedSignal = StatUtil.gaussianSmoother(composite, smoothingWidth);
				double[] smoothedControl = StatUtil.gaussianSmoother(contComposite, smoothingWidth);

				double ratio = getMaxCounts(smoothedSignal) * getMinCounts(smoothedControl)
						/ (getMinCounts(smoothedSignal) * getMaxCounts(smoothedControl));

				writer.println("Sample Name " + rep.getSignal().getName() + "\tRatio: " + ratio + "\tSignal: "
						+ rep.getSignal().getHitCount());
			}
		}
		writer.close();
	}

	public double getMaxCounts(double[] array) {
		double maxCounts = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i] > maxCounts)
				maxCounts = array[i];
		}
		return maxCounts;
	}

	public double getMinCounts(double[] array) {
		double minCounts = Double.MAX_VALUE;
		for (int i = 0; i < array.length; i++) {
			if (array[i] < minCounts)
				minCounts = array[i];
		}
		return minCounts;
	}

	public static void main(String[] args) throws FileNotFoundException {
		ArgParser ap = new ArgParser(args);
		if ((!ap.hasKey("peaks") && !ap.hasKey("regions"))) {
			System.err.println("please input peak files and region files.");
			System.err.println("Usage:\n " + "NucleosomeEnrichmentProfiler\n " + "--geninfo <genome info file> \n "
					+ "--expt <file name> AND --ctrl <file name> AND --format <SAM/BAM/BED/IDX>\n "
					+ "--peaks <file containing coordinates of peaks> \n " + "\nOPTIONS:\n "
					+ "--out <output directory (default = working directory)> \n "
					+ "--win <window of reads to take around peaks (default=200)> \n "
					+ "--bai <file path to bai file> \n "
					+ "--readshift <number of base pair for read shift (default=7)>\n "
					+ "--smooth <window of gaussian kernel smoothing (default=10)> \n " + "");
			System.exit(0);
		}

		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);
		econf.setPerBaseReadFiltering(false);
		ExperimentManager manager = new ExperimentManager(econf);

		// parsing command line arguments
		int win = Args.parseInteger(args, "win", 500);
		int fivePrimeShift = Args.parseInteger(args, "readshift", 6);
		int smooth = Args.parseInteger(args, "smooth", 10);
		List<StrandedPoint> spoints = RegionFileUtilities.loadStrandedPointsFromFile(gconf.getGenome(),
				ap.getKeyValue("peaks"));
		// Get outdir and outbase and make them;
		String outbase = Args.parseString(args, "out", System.getProperty("user.dir"));

		FeatureCountsLoader fcLoader = new FeatureCountsLoader(gconf, spoints, win);
		fcLoader.setFivePrimeShift(fivePrimeShift);

		NucleosomeEnrichmentProfiler profile = new NucleosomeEnrichmentProfiler(outbase, econf, manager, fcLoader,
				smooth);
		profile.execute();

		manager.close();
	}
}
