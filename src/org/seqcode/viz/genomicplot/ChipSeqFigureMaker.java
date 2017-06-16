package org.seqcode.viz.genomicplot;

import java.io.File;

import javax.swing.JFrame;

import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.Species;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.paintable.PaintableFrame;

public class ChipSeqFigureMaker {

	private GenomeConfig gconfig = null;
	private ExptConfig econfig = null;
	private FigureOptions options;
	private PaintableFrame plotter;

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if (!ap.hasKey("options")) {
			System.err.println(
					"Usage:\n " + "ChipSeqFigureMaker " + "--options <file name> " + "--species <species;genome> ");
			return;
		}

		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		String ofile = ap.getKeyValue("options");
		ChipSeqFigureMaker figure = new ChipSeqFigureMaker(ofile, gcon, econ);

	}

	public ChipSeqFigureMaker(String optionFile, GenomeConfig g, ExptConfig e) {
		gconfig = g;
		econfig = e;
		options = new FigureOptions(gconfig, econfig);
		options.loadOptions(new File(optionFile));

		// Paint the picture
		ChipSeqFigurePaintable painter = new ChipSeqFigurePaintable(options);
		plotter = new PaintableFrame("Genomic Data", painter);
		plotter.setSize(options.screenSizeX, options.screenSizeY);
		plotter.setVisible(true);
		plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

}
