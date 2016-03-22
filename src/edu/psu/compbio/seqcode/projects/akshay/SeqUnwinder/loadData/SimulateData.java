package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.loadData;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;
import edu.psu.compbio.seqcode.projects.akshay.utils.SimulateBindingSite;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

public class SimulateData {
	/** Total number of training examples */
	public double N = 10000.0;
	public int[] subclasses = {1700,1200,1100,1400,1200,3400};
	public String[] classNames = {"I;B","I;A","I;C","II;C","II;A","II;B"};
	/** Motifs to inset */
	public Map<String,WeightMatrix> motifs = new HashMap<String,WeightMatrix>();
	public GenomeConfig gcon;
	public MarkovBackgroundModel markov;
	public int win=150;
	
	public SimulateData(GenomeConfig g) {
		gcon = g;
	}
	
	//Settors
	public void setMotifs(Map<String,WeightMatrix> ms){motifs = ms;}
	public void setN(double n){N=n;}
	public void setBack(MarkovBackgroundModel back){markov = back;}
	
	public void execute(){
		// Generate simulated binding sites for each class
		StringBuilder output = new StringBuilder();
		for(int i=0; i<classNames.length; i++){
			String[] labels = classNames[i].split(";");
			// Add all the motifs for this class into a list
			SimulateBindingSite bs = new SimulateBindingSite(gcon);
			
			List<WeightMatrix> mots = new ArrayList<WeightMatrix>();
			mots.add(motifs.get(labels[0]));mots.add(motifs.get(labels[1]));
			List<Double> insrtRate = new ArrayList<Double>();
			insrtRate.add(1.0);insrtRate.add(1.0);
			bs.setBack(markov);
			bs.setInsertRate(insrtRate);
			bs.setLen(win);
			bs.setMotifs(mots);
			
			List<String> seq = bs.execute();
			for(String s : seq){
				output.append(s);output.append("\t");output.append(classNames[i]);output.append("\n");
			}
			
		}
		
		output.deleteCharAt(output.length()-1);
		System.out.println(output.toString());
		
	}
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws IOException, ParseException{
		GenomeConfig gc = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		
		// Set motifs
		String motifFilename = ap.getKeyValue("motiffile");
		if(!ap.hasKey("motiffile")){
			System.err.println("Please provide motifs that need to be inserted!!");
			System.exit(0);
		}
		
		String backFile =ap.getKeyValue("back");
		if(!ap.hasKey("back")){
			System.err.println("Please provide a background model file!!");
			System.exit(0);
		}
		
		MarkovBackgroundModel back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gc.getGenome());
		
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifImport.setBackground(back);
		List<WeightMatrix> mots = new ArrayList<WeightMatrix>();
		mots.addAll(motifImport.readTransfacMatrices(motifFilename));
		Map<String,WeightMatrix> motsWithLabs = new HashMap<String,WeightMatrix>();
		motsWithLabs.put("I", mots.get(0));
		motsWithLabs.put("II", mots.get(1));
		motsWithLabs.put("A", mots.get(2));
		motsWithLabs.put("B", mots.get(3));
		motsWithLabs.put("C", mots.get(4));
		
		SimulateData sd = new SimulateData(gc);
		sd.setMotifs(motsWithLabs);
		sd.setBack(back);
		sd.execute();
		
	}
	
	
}
