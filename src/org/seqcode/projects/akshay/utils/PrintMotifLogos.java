package org.seqcode.projects.akshay.utils;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.genome.GenomeConfig;
import org.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import org.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.BackgroundModelIO;
import org.seqcode.projects.multigps.utilities.Utils;
import org.seqcode.projects.shaun.FreqMatrixImport;


public class PrintMotifLogos {
	
	protected List<WeightMatrix> motifs=new ArrayList<WeightMatrix>();
	
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
