package org.seqcode.deepseq.events;

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

public class BindingModelMaker {
	protected ExperimentManager experiments;
	protected List<StrandedPoint> points;
	protected Integer win;

	public BindingModelMaker(List<StrandedPoint> p, ExperimentManager e, int w) {
		experiments = e;
		win = w;
		points = p;
	}

	public void execute() {
		double[] watson = new double[win];
		double[] crick = new double[win];

		for (ExperimentCondition c : experiments.getConditions()) {
			for (ControlledExperiment r : c.getReplicates()) {
				// Reset
				for (int w = 0; w < win; w++) {
					watson[w] = 0;
					crick[w] = 0;
				}

				// Iterate through points
				for (StrandedPoint pt : points) {
					// Load reads
					List<StrandedBaseCount> wReads = r.getSignal().getStrandedBases(pt.expand(win), pt.getStrand());
					List<StrandedBaseCount> cReads = r.getSignal().getStrandedBases(pt.expand(win),
							pt.getStrand() == '+' ? '-' : '+');

					if (pt.getStrand() == '+') {
						for (StrandedBaseCount sbc : wReads) {
							int sdist = sbc.getCoordinate() - pt.getLocation() + (win / 2);
							if (sdist >= 0 && sdist < win)
								watson[sdist] += sbc.getCount();
						}
						for (StrandedBaseCount sbc : cReads) {
							int sdist = pt.getLocation() - sbc.getCoordinate() + (win / 2);
							if (sdist >= 0 && sdist < win)
								crick[sdist] += sbc.getCount();
						}
					} else {
						for (StrandedBaseCount sbc : wReads) {
							int sdist = pt.getLocation() - sbc.getCoordinate() + (win / 2);
							if (sdist >= 0 && sdist < win)
								watson[sdist] += sbc.getCount();
						}
						for (StrandedBaseCount sbc : cReads) {
							int sdist = sbc.getCoordinate() - pt.getLocation() + (win / 2);
							if (sdist >= 0 && sdist < win)
								crick[sdist] += sbc.getCount();
						}
					}
				}

				// Print
				for (int w = 0; w < win; w++) {
					System.out.println(w - (win / 2) + "\t" + watson[w] + "\t" + crick[w]);
				}
			}
		}
	}

	// Main
	public static void main(String[] args) {
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if (gcon.helpWanted()) {
			System.err.println("BindingModelMaker:");
			System.err.println("\t--points <stranded point file>");
			System.err.println("\t--win <window around points>");
		} else {
			ExperimentManager manager = new ExperimentManager(econ);

			int w = Args.parseInteger(args, "win", 400);
			String pFile = Args.parseString(args, "points", null);
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);
			BindingModelMaker maker = new BindingModelMaker(pts, manager, w);
			maker.execute();

			manager.close();
		}
	}
}
