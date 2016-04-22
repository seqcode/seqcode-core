package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

public class PrintMotifLogos {
	
	protected List<WeightMatrix> motifs;
	
	// Load freq matrices
	public void loadMotifsFromFile(String filename, MarkovBackgroundModel b) {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifImport.setBackground(b);
		motifs.addAll(motifImport.readTransfacMatrices(filename));
	}
	
	public void execute(){
		// Finally, draw the motif logos
		for(WeightMatrix fm : motifs){
			File motifFileName = new File(fm.getName()+".png");
			Utils.printMotifLogo(fm, motifFileName, 150, fm.getName(), true);
		}
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		
		ArgParser ap = new ArgParser(args);
		GenomeConfig gcon = new GenomeConfig(args);
		String motfile = ap.getKeyValue("motfile");
		String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;

		MarkovBackgroundModel back;

		if(backFile == null){
			back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gcon.getGenome()));
		}else{
			back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gcon.getGenome());
		}
		
		PrintMotifLogos plotter = new PrintMotifLogos();
		plotter.loadMotifsFromFile(motfile, back);
		
		plotter.execute();

	}
	
	

}
