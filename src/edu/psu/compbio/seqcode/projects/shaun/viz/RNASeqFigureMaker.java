package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.awt.Color;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableFrame;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintablePanel;

public class RNASeqFigureMaker {

	private Genome gen=null;
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
        
        try {
        	Pair<Species, Genome> pair = Args.parseGenome(args);
        	Species currorg = pair.car();
        	Genome currgen = pair.cdr();
            String ofile = ap.getKeyValue("options");
            RNASeqFigureMaker figure = new RNASeqFigureMaker(ofile, currgen);
            
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public RNASeqFigureMaker(String optionFile, Genome g){
		gen=g;
		options = new FigureOptions(g);
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
