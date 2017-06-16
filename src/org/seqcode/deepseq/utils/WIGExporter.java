package org.seqcode.deepseq.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import org.seqcode.deepseq.ReadHit;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;

/**
 * Outputs a fixed-step WIG format file for a deep-seq experiment.
 * 
 * @author Shaun Mahony
 * @version %I%, %G%
 */
public class WIGExporter {
	private Species org;
	private Genome gen;
	private Sample sample;
	protected int[] stackedHitCounts;
	private int winSize = 20, winStep = 20;
	private int readLength = 1, read5PrimeExt = 0, read3PrimeExt = 200;
	private String outName = "out";
	private String trackName = "out";
	private String trackDesc = "out";
	private String trackColor = "0,0,255";
	private int trackYMax = -1;
	private int perBaseMax = -1;
	private boolean needlefiltering = false;
	private boolean cacheAllData;

	public static void main(String[] args) throws SQLException, NotFoundException {

		WIGExporter wig = new WIGExporter(args);
		wig.execute();
		wig.close();

		System.exit(0);
	}

	public WIGExporter(String[] args) {
		if (args.length == 0) {
			System.err.println("WIGExporter usage:\n" + "\t--species <organism;genome>\n"
					+ "\t--(rdb)expt <experiment names>\n" + "\t--read5ext <Extension on the 5' end>\n"
					+ "\t--read3ext <Extension on the 3' end>\n" + "\t--pbmax <max read count per base>\n"
					+ "\t--winsize <window size/step in WIG file>\n" + "\t--name <string to use as track name>\n"
					+ "\t--description <string to use as track description>\n" + "\t--ylimit <default track y max>\n"
					+ "\t--color <R,G,B>\n" + "\t--out <output file name>");
			System.exit(1);
		}
		ArgParser ap = new ArgParser(args);
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig c = new ExptConfig(gcon.getGenome(), args);

		// WIGExporter specific options
		outName = Args.parseString(args, "out", outName);
		trackName = Args.parseString(args, "name", trackName);
		trackDesc = Args.parseString(args, "description", trackDesc);
		trackColor = Args.parseString(args, "color", trackColor);
		read5PrimeExt = Args.parseInteger(args, "read5ext", read5PrimeExt);
		read3PrimeExt = Args.parseInteger(args, "read3ext", read3PrimeExt);
		readLength = Args.parseInteger(args, "readlen", readLength);
		winSize = Args.parseInteger(args, "winsize", winSize);

		// General options processed directly by ExptConfig

		perBaseMax = (int) c.getPerBaseMax();
		if (perBaseMax != -1) {
			needlefiltering = true;
		}
		if (ap.hasKey("ylimit")) {
			trackYMax = Args.parseInteger(args, "ylimit", -1);
		}
		winStep = winSize;
		ExperimentManager manager = new ExperimentManager(c);
		sample = manager.getSamples().get(0);
		cacheAllData = c.getCacheAllData();
		gen = gcon.getGenome();

		// Load the experiments
		// List<SeqLocator> dbexpts = Args.parseSeqExpt(args, "dbexpt");
		// List<SeqLocator> rdbexpts = Args.parseSeqExpt(args,"rdbexpt");
		// List<File> expts = Args.parseFileHandles(args, "expt");
		// boolean nonUnique = ap.hasKey("nonunique") ? true : false;
		// String fileFormat = Args.parseString(args, "format", "ELAND");
		// if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
		// expt= new DeepSeqExpt(gen, expts, nonUnique, fileFormat,
		// (int)readLength);
		// }else if (dbexpts.size() > 0 && expts.size() == 0) {
		// expt = new DeepSeqExpt(gen, dbexpts, "db", (int)readLength);
		// dbconnected = true;
		// }else if (rdbexpts.size()>0 && expts.size() == 0){
		// expt = new DeepSeqExpt(gen, rdbexpts, "readdb", -1);
		// dbconnected=true;
		// }else {
		// logger.error("Must provide either an aligner output file or Gifford
		// lab DB experiment name for the signal experiment (but not both)");
		// System.exit(1);
		// }
		// logger.info("Expt hit count: " + (int) expt.getHitCount() + ",
		// weight: " + (int) expt.getWeightTotal());

		read3PrimeExt = Math.max(0, read3PrimeExt - readLength);

		sample.initializeCache(cacheAllData, null);
	}

	public void close() {
		sample.close();
	}

	public void execute() {
		try {
			FileWriter fw = new FileWriter(outName + ".wig");

			double basesDone = 0, printStep = 10000000, numPrint = 0;
			if (trackName.equals("out"))
				trackName = outName;
			if (trackDesc.equals("out"))
				trackDesc = outName;

			// Print the header
			fw.write("track type=wiggle_0 name=\"" + trackName + "\" description=\"" + trackDesc + " summary\""
					+ " visibility=full color=" + trackColor + " ");
			if (trackYMax > 0)
				fw.write("autoScale=off viewLimits=0:" + trackYMax + " ");
			fw.write("\n");

			ChromRegionIterator chroms = new ChromRegionIterator(gen);
			while (chroms.hasNext()) {
				NamedRegion currentRegion = chroms.next();

				// Split the job up into chunks of 100Mbp
				for (int x = currentRegion.getStart(); x <= currentRegion.getEnd(); x += 100000000) {
					int y = x + 100000000;
					if (y > currentRegion.getEnd()) {
						y = currentRegion.getEnd();
					}
					Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);

					ArrayList<ReadHit> hits = new ArrayList<ReadHit>();
					hits.addAll(sample.exportExtReadHits(currSubRegion, readLength, 0, read5PrimeExt, read3PrimeExt));
					double stackedHitCounts[] = makeHitLandscape(hits, currSubRegion, perBaseMax, '.');

					boolean recording = false;
					// Scan regions
					for (int i = currSubRegion.getStart(); i < currSubRegion.getEnd()
							- (int) winSize; i += (int) winStep) {
						Region currWin = new Region(gen, currentRegion.getChrom(), i, (int) (i + winSize - 1));

						int binid = (int) Math.max(0,
								((double) (currWin.getStart() - currSubRegion.getStart()) / winStep));
						double winHits = (double) stackedHitCounts[binid];

						if (winHits > 0) {
							if (!recording) {
								fw.write("fixedStep chrom=chr" + currSubRegion.getChrom() + " start=" + (i + 1)
										+ " step=" + winStep + " span=" + winSize + "\n");
								recording = true;
							}
							fw.write(String.format("%.1f\n", winHits));
						} else {
							recording = false;
						}

						// Print out progress
						basesDone += winStep;
						if (basesDone > numPrint * printStep) {
							if (numPrint % 10 == 0) {
								System.out.print(String.format("(%.0f)", (numPrint * printStep)));
							} else {
								System.out.print(".");
							}
							if (numPrint % 50 == 0 && numPrint != 0) {
								System.out.print("\n");
							}
							numPrint++;
						}
					}
				}
			}
			System.out.print("\n");
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected double[] makeHitLandscape(ArrayList<ReadHit> hits, Region currReg, int perBaseMax, char strand) {
		int numBins = (int) (currReg.getWidth() / winStep);
		double[] landscape = new double[numBins + 1];
		double[] startcounts = new double[(int) currReg.getWidth() + 1];
		for (int i = 0; i <= numBins; i++) {
			landscape[i] = 0;
		}
		for (int i = 0; i <= currReg.getWidth(); i++) {
			startcounts[i] = 0;
		}
		for (ReadHit r : hits) {
			if (strand == '.' || r.getStrand() == strand) {
				int offset5 = inBounds(r.getFivePrime() - currReg.getStart(), 0, currReg.getWidth());
				if (!needlefiltering || (startcounts[offset5] <= perBaseMax)) {
					int binstart = inBounds((int) ((double) (r.getStart() - currReg.getStart()) / winStep), 0, numBins);
					int binend = inBounds((int) ((double) (r.getEnd() - currReg.getStart()) / winStep), 0, numBins);
					for (int i = binstart; i <= binend; i++) {
						landscape[i] += r.getWeight();
					}
					if (needlefiltering && (startcounts[offset5] + r.getWeight() > perBaseMax))
						startcounts[offset5] = perBaseMax;
					else
						startcounts[offset5] += r.getWeight();
				}
			}
		}
		return (landscape);
	}

	// keep the number in bounds
	protected final double inBounds(double x, double min, double max) {
		if (x < min) {
			return min;
		}
		if (x > max) {
			return max;
		}
		return x;
	}

	protected final int inBounds(int x, int min, int max) {
		if (x < min) {
			return min;
		}
		if (x > max) {
			return max;
		}
		return x;
	}
}
