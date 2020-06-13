package org.seqcode.viz.genomicplot;

import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;

import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.viz.paintable.PaintableFrame;


public class ChipSeqFigureMaker {

	private GenomeConfig gconfig=null;
	private ExptConfig econfig=null;
	private FigureOptions options; 
	private PaintableFrame plotter;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("options")) { 
            System.err.println("Usage:\n " +
                               "ChipSeqFigureMaker\n" +
                               "\t--options <file name>\n" +
                               "\t--species <species;genome>\t"+
                               "\t--batch [batch mode]\n" +
                               "\t--out <file name for batch mode>");
            return;
        }
        
        GenomeConfig gcon = new GenomeConfig(args);
    	ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
        String ofile = ap.getKeyValue("options");
        File outfile = new File(ap.getKeyValue("out"));
        boolean isBatch = ap.hasKey("batch");
        ChipSeqFigureMaker figure = new ChipSeqFigureMaker(ofile, gcon, econ, isBatch, outfile);
        
	}
	
	public ChipSeqFigureMaker(String optionFile, GenomeConfig g, ExptConfig e, boolean isBatch, File outFile){
		gconfig=g;
		econfig=e;
		options = new FigureOptions(gconfig, econfig);
		options.loadOptions(new File(optionFile));
		//Paint the picture
		ChipSeqFigurePaintable painter = new ChipSeqFigurePaintable(options);
		
		if(isBatch){
			System.setProperty("java.awt.headless", "true");
			try {
				painter.saveImage(outFile, options.screenSizeX, options.screenSizeY, true);
				System.exit(0);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}else{
			plotter = new PaintableFrame("Genomic Data", painter);
			plotter.setSize(options.screenSizeX, options.screenSizeY);
			plotter.setVisible(true);
			plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);		
		}
	}

}
