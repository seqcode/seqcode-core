package org.seqcode.deepseq.experiments;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.data.seqdata.SeqDataLoader;
import org.seqcode.deepseq.hitloaders.HitLoader;
import org.seqcode.deepseq.hitloaders.HitLoaderFactory;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.Pair;

/**
 * SampleLoader allows ad hoc loading of a set of individual Samples and
 * associated HitLoaders. No organization/relationship is assumed or recorded
 * between the loaded Samples. ExperimentManager is the preferred way to load a
 * set of experiments. In rare cases, an application may just need to load one
 * or two Samples without requiring the overhead that comes with the Experiment
 * tree. That's where this class comes in.
 * 
 * @author mahony
 */
public class SampleLoader {
	protected ExptConfig econfig;
	protected HitLoaderFactory hlfactory;
	protected SeqDataLoader sdloader = null;
	protected Genome gen;
	protected boolean loadReads = true;

	// Experiment tree elements
	protected HashMap<String, HitLoader> loaders = new HashMap<String, HitLoader>();
	protected List<Sample> samples = new ArrayList<Sample>();
	protected Map<String, Sample> samplesByName = new HashMap<String, Sample>();

	/**
	 * Constructor: Using arguments loaded by the ExptConfig, initialize (in
	 * this order): HitLoaders, Samples.
	 * 
	 * @param c
	 *            : ExptConfig
	 * @param loadReads
	 *            : boolean. for some applications, reads do not have to be
	 *            loaded. Use with caution.
	 */
	public SampleLoader(ExptConfig c) {
		this(c, true);
	}

	public SampleLoader(ExptConfig c, boolean loadReads) {
		econfig = c;
		gen = econfig.getGenome();
		hlfactory = new HitLoaderFactory(econfig);

		this.loadReads = loadReads;

		List<ExptDescriptor> descriptors = econfig.getExperimentDescriptors();
		if (descriptors != null && descriptors.size() > 0)
			loadSamples(descriptors);
	}

	/**
	 * Load some samples
	 * 
	 * @param descriptors
	 */
	public void loadSamples(List<ExptDescriptor> descriptors) {
		HashMap<String, Sample> allSamples = new HashMap<String, Sample>();
		int sampCount = 0;

		// Pre-step; do we need a SeqDataLoader?
		boolean makeSeqDataLoader = false;
		for (ExptDescriptor e : descriptors) {
			for (Pair<String, String> source : e.sources) {
				String type = source.cdr();
				if (type.equals("READDB"))
					makeSeqDataLoader = true;
			}
		}
		if (makeSeqDataLoader)
			try {
				sdloader = new SeqDataLoader();
			} catch (SQLException e1) {
				e1.printStackTrace();
			} catch (IOException e1) {
				e1.printStackTrace();
			}

		// Firstly, initialize all hit loaders.
		// This is done in a separate first pass, because it is possible (albeit
		// unlikely)
		// that multiple conditions share the same hit loader, and you don't
		// want to load things twice.
		for (ExptDescriptor e : descriptors) {
			if (econfig.getPrintLoadingProgress())
				System.err.println("Processing HitLoaders for:\t" + e.condition + "\t" + e.replicate);
			for (Pair<String, String> source : e.sources) {
				String name = source.car();
				String type = source.cdr();
				if (type.equals("READDB")) { // ReadDB HitLoader
					if (!loaders.containsKey(name)) {
						HitLoader hl = hlfactory.makeReadDBHitLoader(sdloader, name);
						// hit loader does not have to be sourced here -- that
						// happens in the samples part below
						loaders.put(name, hl);
					}
				} else { // Assume File HitLoader
					if (!loaders.containsKey(name)) {
						HitLoader hl = hlfactory.makeFileHitLoader(name, type, econfig.getNonUnique());
						// hit loader does not have to be sourced here -- that
						// happens in the samples part below
						loaders.put(name, hl);
					}
				}
			}
		}

		// Secondly, load the samples (load each sample name once)
		for (ExptDescriptor e : descriptors) {
			String sampleName = e.getName();

			if (econfig.getPrintLoadingProgress() && loadReads)
				System.err.print("Loading data from " + sampleName);

			if (!allSamples.containsKey(sampleName)) {
				Sample samp = new Sample(sampCount, econfig, sampleName, e.perBaseMaxReads, e.signal);
				allSamples.put(sampleName, samp);
				samples.add(samp);
				samplesByName.put(sampleName, samp);
				sampCount++;
			}
			for (Pair<String, String> source : e.sources) {
				String name = source.car();
				allSamples.get(sampleName).addHitLoader(loaders.get(name));
			}
			if (loadReads) {
				allSamples.get(sampleName).initializeCache(econfig.getCacheAllData(),
						econfig.getInitialCachedRegions());
				if (econfig.getPrintLoadingProgress())
					System.err.println(String.format("\tLoaded:\t%.1f", allSamples.get(sampleName).getHitCount()));
			}
		}
		// Merge estimated genomes if necessary (v. messy if Samples are loaded
		// during an analysis...)
		if (gen == null) {
			List<Genome> estGenomes = new ArrayList<Genome>();
			for (String s : allSamples.keySet())
				estGenomes.add(allSamples.get(s).getGenome());
			gen = econfig.mergeEstGenomes(estGenomes);
			for (String s : allSamples.keySet())
				allSamples.get(s).setGenome(gen);
		}

		if (sdloader != null) {
			sdloader.close();
			sdloader = null;
		}
	}

	// Accessors
	public List<Sample> getSamples() {
		return samples;
	}

	public Sample getSample(String s) {
		return samplesByName.get(s);
	}

	/**
	 * Call any cleanup methods
	 */
	public void close() {
		for (String l : loaders.keySet()) {
			loaders.get(l).cleanup();
		}
		for (Sample s : samples) {
			s.close();
		}
	}

	/**
	 * This main method is only for testing the ExperimentManager system
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		if (econfig.helpWanted()) {
			System.err.println("ExperimentManager debugging:");
			System.err.println(econfig.getArgsList());
		} else {
			ExperimentManager manager = new ExperimentManager(econfig);

			System.err.println("ExptTypes:\t" + manager.getExptTypes().size());
			for (ExperimentType t : manager.getExptTypes()) {
				System.err
						.println("ExptType " + t.getName() + ":\t#Experiments:\t" + t.getExptTypeExperiments().size());
			}
			System.err.println("ExptTargets:\t" + manager.getTargets().size());
			for (ExperimentTarget t : manager.getTargets()) {
				System.err.println("Target " + t.getName() + ":\t#Experiments:\t" + t.getTargetExperiments().size());
			}

			manager.close();
		}
	}

}
