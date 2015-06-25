package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.io.File;

import javax.swing.JFrame;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableFrame;

public class ChipSeqFigureMaker {

	private Genome gen=null;
	private FigureOptions options; 
	private PaintableFrame plotter;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("options")) { 
            System.err.println("Usage:\n " +
                               "ChipSeqFigureMaker " +
                               "--options <file name> " +
                               "--species <species;genome> ");
            return;
        }
        
        try {
        	Pair<Species, Genome> pair = Args.parseGenome(args);
        	Species currorg = pair.car();
        	Genome currgen = pair.cdr();
            String ofile = ap.getKeyValue("options");
            ChipSeqFigureMaker figure = new ChipSeqFigureMaker(ofile, currgen);
            
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public ChipSeqFigureMaker(String optionFile, Genome g){
		gen=g;
		options = new FigureOptions(g);
		options.loadOptions(new File(optionFile));
		
		//Paint the picture
		ChipSeqFigurePaintable painter = new ChipSeqFigurePaintable(options);
		plotter = new PaintableFrame("Genomic Data", painter);
		plotter.setSize(options.screenSizeX, options.screenSizeY);
		plotter.setVisible(true);
		plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);		
	}

}
