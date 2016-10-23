package org.seqcode.viz.genomicplot;

import java.awt.Color;
import java.io.File;

import javax.swing.JFrame;

import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.viz.paintable.PaintableFrame;


public class RNASeqFigureMaker {

	private GenomeConfig gconfig=null;
	private ExptConfig econfig=null;
	private FigureOptions options;
	private PaintableFrame plotter;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("options")) { 
            System.err.println("Usage:\n " +
                               "RNASeqFigureMaker " +
                               "--options <file name> " +
                               "--species <species;genome> ");
            return;
        }
        
    	GenomeConfig gcon = new GenomeConfig(args);
    	ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
        String ofile = ap.getKeyValue("options");
        RNASeqFigureMaker figure = new RNASeqFigureMaker(ofile, gcon, econ);
    }
	
	public RNASeqFigureMaker(String optionFile, GenomeConfig g, ExptConfig e){
		gconfig=g;
		econfig=e;
		options = new FigureOptions(gconfig, econfig);
		options.loadOptions(new File(optionFile));
		
		//Paint the picture
		RNASeqFigurePaintable painter = new RNASeqFigurePaintable(options);
		plotter = new PaintableFrame("Genomic Data", painter);
		plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		plotter.setBackground(Color.white);
		plotter.setSize(options.screenSizeX, options.screenSizeY);
		plotter.setVisible(true);

	}

}
