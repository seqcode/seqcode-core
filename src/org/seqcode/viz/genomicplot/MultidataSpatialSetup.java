package org.seqcode.viz.genomicplot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.swing.JFrame;

import org.seqcode.data.seqdata.SeqLocator;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.ExptDescriptor;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.experiments.SampleLoader;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.ExonicGene;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gsebricks.verbs.location.PointParser;
import org.seqcode.gsebricks.verbs.location.RegionParser;
import org.seqcode.gsebricks.verbs.location.StrandedRegionParser;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.paintable.PaintableFrame;

//Note: RAR figure was defined for mm8
public class MultidataSpatialSetup {

	private MultidataSpatialPaintable painter;
	private PaintableFrame plotter;
	private int screenSizeX = 1000, screenSizeY = 900;
	private String inputFile;
	private Map<String, ExptDescriptor> experiments = new HashMap<String, ExptDescriptor>();
	private List<Sample> times = new ArrayList<Sample>();
	private List<String> timeLabels = new ArrayList<String>();
	private List<ExonicGene> genes = new ArrayList<ExonicGene>();
	private Map<String, List<Double>> expression = new HashMap<String, List<Double>>();
	private Map<String, List<Region>> sites = new HashMap<String, List<Region>>();
	private List<Pair<String, StrandedRegion>> motifs = new ArrayList<Pair<String, StrandedRegion>>();
	private List<Region> lits = new ArrayList<Region>();
	private int rstart = 0;
	private int rstop = 0;
	private String chr;
	public GenomeConfig gconfig;
	public Genome gen;
	public ExptConfig econfig;
	public SampleLoader sampleLoader;

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if (!ap.hasKey("data")) {
			System.err.println("Usage:\n " + "MultidataSpatialSetup " + "--data <file name> ");
			return;
		}
		String dfile = ap.getKeyValue("data");
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);

		MultidataSpatialSetup setup = new MultidataSpatialSetup(dfile, gcon, econ);
	}

	public MultidataSpatialSetup(String df, GenomeConfig gcon, ExptConfig econ) {
		gconfig = gcon;
		econfig = econ;
		sampleLoader = new SampleLoader(econfig);
		inputFile = df;

		// Load the file contents
		loadFile(inputFile);

		Region gRegion = new Region(gen, chr, rstart, rstop);

		// Paint the picture
		MultidataSpatialPaintable painter = new MultidataSpatialPaintable(timeLabels, times, gRegion, genes, expression,
				sites, motifs, lits);
		plotter = new PaintableFrame("Genomic Data", painter);
		plotter.setSize(screenSizeX, screenSizeY);
		plotter.setVisible(true);
		plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	private void loadFile(String inF) {
		try {
			File aFile = new File(inF);
			if (aFile.isFile()) {
				BufferedReader reader;
				reader = new BufferedReader(new FileReader(aFile));
				String line;
				while ((line = reader.readLine()) != null) {
					String[] tokens = line.split("\\t");
					if (tokens[0].equals("timepoints")) {
						for (int i = 1; i < tokens.length; i++) {
							timeLabels.add(tokens[i]);
						}
					} else if (tokens.length == 2) {// Setup annotations
						if (tokens[0].equals("regionstart")) {
							rstart = new Integer(tokens[1]).intValue();
						}
						if (tokens[0].equals("regionstop")) {
							rstop = new Integer(tokens[1]).intValue();
						}
						if (tokens[0].equals("regionchr")) {
							chr = new String(tokens[1]);
						}
					} else if (tokens.length >= 15) {// A gene
						int start = new Integer(tokens[4]).intValue();
						int stop = new Integer(tokens[5]).intValue();
						String name = tokens[12];
						String id = tokens[1];
						char str = tokens[3].charAt(0);
						String[] estarts = tokens[9].split(",");
						String[] estops = tokens[10].split(",");
						ExonicGene e = new ExonicGene(gen, chr, start, stop, name, id, str, "Manual");
						for (int i = 0; i < estarts.length; i++) {
							Region r = new Region(gen, chr, new Integer(estarts[i]).intValue(),
									new Integer(estops[i]).intValue());
							e.addExon(r);
						}
						genes.add(e);
						System.out.println(name);
					} else if (tokens[0].equals("Expr")) {// An expression set
						String name = tokens[1];
						ArrayList<Double> vals = new ArrayList<Double>();
						for (int i = 2; i < tokens.length; i++) {
							vals.add(new Double(tokens[i]));
						}
						expression.put(name, vals);
					} else if (tokens[0].equals("Site")) {// A binding site
															// (timepoint,
															// position, zscore)
						String time = tokens[1];
						RegionParser rparser = new RegionParser(gen);
						Region r = rparser.execute(tokens[2]);
						if (!sites.containsKey(time)) {
							sites.put(time, new ArrayList<Region>());
						}
						sites.get(time).add(r);
					} else if (tokens[0].equals("Motif")) {// Motif hits
						StrandedRegionParser parser = new StrandedRegionParser(gen);
						StrandedRegion s = parser.execute(tokens[2]);
						String type = tokens[1];
						if (s != null)
							motifs.add(new Pair<String, StrandedRegion>(type, s));
					} else if (tokens[0].equals("Lit")) {// Known sites
						Region site = new Region(gen, tokens[3], new Integer(tokens[4]).intValue(),
								new Integer(tokens[5]).intValue());
						lits.add(site);
					} else if (tokens[0].equals("Expt")) {// Experiments
						if (!experiments.containsKey(tokens[1])) {
							System.err.println("Experiment: " + tokens[1]);
							experiments.put(tokens[1],
									new ExptDescriptor("", "", tokens[1], "A", true,
											new Pair<String, String>(tokens[2] + ";" + tokens[3], "READDB"),
											econfig.getPerBaseMax()));
						}
					}

				}
				reader.close();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
