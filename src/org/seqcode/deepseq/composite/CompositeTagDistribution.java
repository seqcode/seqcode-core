package org.seqcode.deepseq.composite;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;

/**
 * CompositeTagDistribution: watson/crick tag distributions for a collection of
 * aligned points and the resulting composite
 *
 * Coordinates are 0-based, and the center of the distribution/alignment is
 * defined by the centerOffset variable
 * 
 * @author mahony
 *
 */
public class CompositeTagDistribution {
	protected ExperimentManager exptMan;
	protected List<StrandedPoint> points;
	protected int win;
	protected int centerOffset;
	protected int numConditions;
	protected int numPoints;
	protected double[][] watson; // per-condition watson tags {condition,
									// location}
	protected double[][] crick; // per-condition crick tags {condition,
								// location}
	protected double[][][] perPointWatson; // per-point, per-condition watson
											// tags {point, condition, location}
	protected double[][][] perPointCrick; // per-point, per-condition crick tags
											// {point, condition, location}
	protected HashMap<StrandedPoint, Integer> pointIndex = new HashMap<StrandedPoint, Integer>();
	protected boolean isSignal;

	public CompositeTagDistribution(List<StrandedPoint> points, ExperimentManager eMan, int win, boolean loadSignal) {
		exptMan = eMan;
		this.win = win;
		centerOffset = win / 2;
		this.numConditions = exptMan.getNumConditions();
		this.points = points;
		numPoints = points.size();
		isSignal = loadSignal;

		watson = new double[numConditions][win];
		crick = new double[numConditions][win];
		perPointWatson = new double[numPoints][numConditions][win];
		perPointCrick = new double[numPoints][numConditions][win];

		for (int p = 0; p < numPoints; p++)
			pointIndex.put(points.get(p), p);

		// Reset
		for (int c = 0; c < numConditions; c++) {
			for (int w = 0; w < win; w++) {
				watson[c][w] = 0;
				crick[c][w] = 0;
			}
			for (int p = 0; p < numPoints; p++)
				for (int w = 0; w < win; w++) {
					perPointWatson[p][c][w] = 0;
					perPointCrick[p][c][w] = 0;
				}
		}

		for (ExperimentCondition cond : exptMan.getConditions()) {
			for (ControlledExperiment rep : cond.getReplicates()) {

				if (loadSignal || rep.hasControl()) {
					// Iterate through points
					int p = 0;
					for (StrandedPoint pt : points) {
						// Load reads
						List<StrandedBaseCount> wReads = loadSignal
								? rep.getSignal().getStrandedBases(pt.expand(win), pt.getStrand())
								: rep.getControl().getStrandedBases(pt.expand(win), pt.getStrand());
						List<StrandedBaseCount> cReads = loadSignal
								? rep.getSignal().getStrandedBases(pt.expand(win), pt.getStrand() == '+' ? '-' : '+')
								: rep.getControl().getStrandedBases(pt.expand(win), pt.getStrand() == '+' ? '-' : '+');

						if (pt.getStrand() == '+') {
							for (StrandedBaseCount sbc : wReads) {
								int sdist = sbc.getCoordinate() - pt.getLocation() + (win / 2);
								if (sdist >= 0 && sdist < win) {
									watson[cond.getIndex()][sdist] += sbc.getCount();
									perPointWatson[p][cond.getIndex()][sdist] += sbc.getCount();
								}
							}
							for (StrandedBaseCount sbc : cReads) {
								int sdist = sbc.getCoordinate() - pt.getLocation() + (win / 2);
								if (sdist >= 0 && sdist < win) {
									crick[cond.getIndex()][sdist] += sbc.getCount();
									perPointCrick[p][cond.getIndex()][sdist] += sbc.getCount();
								}
							}
						} else {
							for (StrandedBaseCount sbc : wReads) {
								int sdist = pt.getLocation() - sbc.getCoordinate() + (win / 2);
								if (sdist >= 0 && sdist < win) {
									watson[cond.getIndex()][sdist] += sbc.getCount();
									perPointWatson[p][cond.getIndex()][sdist] += sbc.getCount();
								}
							}
							for (StrandedBaseCount sbc : cReads) {
								int sdist = pt.getLocation() - sbc.getCoordinate() + (win / 2);
								if (sdist >= 0 && sdist < win) {
									crick[cond.getIndex()][sdist] += sbc.getCount();
									perPointCrick[p][cond.getIndex()][sdist] += sbc.getCount();
								}
							}
						}

						p++;
					}
				}
			}
			// Normalize
			double wsum = 0, csum = 0;
			for (int w = 0; w < win; w++) {
				wsum += watson[cond.getIndex()][w];
				csum += crick[cond.getIndex()][w];
			}
			for (int w = 0; w < win; w++) {
				watson[cond.getIndex()][w] /= wsum;
				crick[cond.getIndex()][w] /= csum;
			}
		}
	}

	// Accessors
	public int getWinSize() {
		return win;
	}

	public int getCenterOffset() {
		return centerOffset;
	}

	public int getNumConditions() {
		return numConditions;
	}

	public double[][] getCompositeWatson() {
		return watson;
	}

	public double[][] getCompositeCrick() {
		return crick;
	}

	public double[] getCompositeWatson(ExperimentCondition c) {
		return watson[c.getIndex()];
	}

	public double[] getCompositeCrick(ExperimentCondition c) {
		return crick[c.getIndex()];
	}

	public double[] getPointWatson(StrandedPoint p, ExperimentCondition c) {
		return perPointWatson[pointIndex.get(p)][c.getIndex()];
	}

	public double[] getPointCrick(StrandedPoint p, ExperimentCondition c) {
		return perPointCrick[pointIndex.get(p)][c.getIndex()];
	}

	public double[][] getPointWatsons(int index) {
		return perPointWatson[index];
	}

	public double[][] getPointCricks(int index) {
		return perPointCrick[index];
	}

	public List<StrandedPoint> getPoints() {
		return points;
	}

	public StrandedPoint getPoint(int i) {
		return points.get(i);
	}

	/**
	 * Per-condition sum of tags in composites
	 * 
	 * @return
	 */
	public double[] getCompositeSums() {
		double[] sums = new double[numConditions];
		for (int c = 0; c < numConditions; c++) {
			for (int i = 0; i < win; i++)
				sums[c] = 0;
			for (int i = 0; i < win; i++) {
				sums[c] += watson[c][i] + crick[c][i];
			}
		}
		return sums;
	}

	public String toString(ExperimentCondition cond) {
		String out = "";
		for (int w = 0; w < win; w++) {
			int pos = (w - centerOffset);
			out = out + pos + "\t" + watson[cond.getIndex()][w] + "\t" + crick[cond.getIndex()][w] + "\n";
		}
		return out;
	}

	// Print probs to a file
	public void printProbsToFile(ExperimentCondition cond, String filename) {
		try {
			FileWriter fout = new FileWriter(filename);
			for (int w = 0; w < win; w++) {
				int pos = (w - centerOffset);
				fout.write(pos + "\t" + watson[cond.getIndex()][w] + "\t" + crick[cond.getIndex()][w] + "\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// Main method to make new composite distributions
	public static void main(String[] args) {
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if (args.length == 0) {
			System.err.println("CompositeTagDistribution:" + "\t--cpoints <stranded point file>"
					+ "\t--cwin <window around points>" + "Genome:" + "\t--species <Species;Genome>\n" + "\tOR\n"
					+ "\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" + "Experiment Design File:\n"
					+ "\t--design <file name>\n");
		} else {
			ExperimentManager manager = new ExperimentManager(econ);

			int w = Args.parseInteger(args, "cwin", 400);
			String pFile = Args.parseString(args, "cpoints", null);
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);

			CompositeTagDistribution maker = new CompositeTagDistribution(pts, manager, w, true);

			for (ExperimentCondition cond : manager.getConditions()) {
				String compositeFileName = "out_composite." + cond.getName() + ".txt";
				maker.printProbsToFile(cond, compositeFileName);
			}
			manager.close();
		}
	}

}
