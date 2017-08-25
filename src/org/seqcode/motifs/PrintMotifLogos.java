package org.seqcode.motifs;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.motifdb.CountsBackgroundModel;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.motifs.FreqMatrixImport;

/**
 * For each input motif, print motif logos.
 */

public class PrintMotifLogos {
	
	protected List<WeightMatrix> motifs=new ArrayList<WeightMatrix>();
	
	// Load freq matrices
	public void loadMotifsFromFile(String filename, MarkovBackgroundModel b) {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifImport.setBackground(b);
		motifs.addAll(motifImport.readTransfacMatrices(filename));
	}
	
	public void execute(boolean drawaxis){
		// Finally, draw the motif logos
		for(WeightMatrix fm : motifs){
			File motifFileName = new File(fm.getName()+".png");
			org.seqcode.motifs.DrawMotifs.printMotifLogo(fm, motifFileName, 150, fm.getName(), drawaxis);
			motifFileName = new File(fm.getName()+"_rc.png");
			org.seqcode.motifs.DrawMotifs.printMotifLogo(WeightMatrix.reverseComplement(fm), motifFileName, 150, fm.getName(), drawaxis);
		}
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		
		ArgParser ap = new ArgParser(args);		
        if(!ap.hasKey("motfile")) { 
        	System.err.println("please input motfile.");
            System.err.println("Usage:\n " +
                               "PrintMotifLogos\n " +
                               "--geninfo <genome info file> \n " +
                               "--expt <file name> AND --ctrl <file name> AND --format <SAM/BAM/BED/IDX>\n " +
                               "--motfile <file containing motifs> \n " +
                               "");
            System.exit(0);
        }

		GenomeConfig gcon = new GenomeConfig(args);
		String motfile = ap.getKeyValue("motfile");
		String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;
		boolean showaxis = ap.hasKey("showaxix");

		MarkovBackgroundModel back;

		if(backFile == null){
			back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gcon.getGenome()));
		}else{
			back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gcon.getGenome());
		}
		
		PrintMotifLogos plotter = new PrintMotifLogos();
		plotter.loadMotifsFromFile(motfile, back);
		
		plotter.execute(showaxis);

	}
	
	

}

