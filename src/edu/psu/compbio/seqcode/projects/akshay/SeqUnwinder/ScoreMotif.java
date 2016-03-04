package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

public class ScoreMotif {
	
	public List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	public int Kmin = 4;
	public int Kmax = 6;
	public int numK;
	public HashMap<String,double[]> kmerweights = new HashMap<String,double[]>();
	/** Model names */
	protected List<String> kmerModelNames = new ArrayList<String>();
	
	
	
	public ScoreMotif() {
		// TODO Auto-generated constructor stub
	}
	
	//Settors
	public void setKmerMin(int m){Kmin = m;}
	public void setKmerMax(int m){Kmax = m;}
	public void setNumK(){
		numK = 0;
		for(int k=Kmin; k<=Kmax; k++ ){
			numK += (int)Math.pow(4, k);
		}
	}
	// Load freq matrices
	public void loadMotifsFromFile(String filename) {
		motifs.addAll(FreqMatrixImport.readTransfacMatrices(filename));
	}

	// Gettors
	public int getKmerBaseInd(String kmer) {
		int baseInd = 0;
		for (int k = Kmin; k < kmer.length(); k++) {
			baseInd += (int) Math.pow(4, k);
		}
		return baseInd;
	}
	
	
	public void execute(){
		StringBuilder sb = new StringBuilder();
		for(WeightMatrix mot : motifs){
			HashSet<String> motKmers = WeightMatrix.getConsensusKmers(mot, Kmin, Kmax);
			for(String modName : kmerModelNames){
				double weight=0.0;
				for(String s : motKmers){
					int ind = getKmerBaseInd(s) + RegionFileUtilities.seq2int(s);
					weight = weight + kmerweights.get(modName)[ind];
				}
				weight = weight/motKmers.size();
				sb.append(mot.getName()+"\t");sb.append(modName+"\t");sb.append(weight);sb.append("\n");
			}
		}
		System.out.println(sb.toString());
	}
	
	
	public void setKmerWeights(String weightFileName) throws NumberFormatException, IOException{
		BufferedReader reader = new BufferedReader(new FileReader(weightFileName));
        String line;
        boolean header = false;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\t");
            if(words[0].charAt(0) == '#' || words[0].contains("Variable") || words[0].contains("Class")){
            	header = true;
            	for(int i=1; i<words.length; i++){
            		kmerModelNames.add(words[i]);
            		kmerweights.put(words[i], new double[numK]);
            	}
            }else{
            	if(!header){
            		System.err.println("Please provide a header in the K-mer weight file");
            		System.exit(1);
            	}
            	int ind = getKmerBaseInd(words[0]) + RegionFileUtilities.seq2int(words[0]);
            	for(int i = 1; i < words.length; i++ ){
            		kmerweights.get(kmerModelNames.get(i-1))[ind] = Double.parseDouble(words[i]);
            	}
            }
           
        }reader.close();
	}
	
	public static void main(String[] args) throws IOException {
		ArgParser ap = new ArgParser(args);
		ScoreMotif runner = new ScoreMotif();

		// Length of the smallest K-mer in the K-mer models
		int minK = Args.parseInteger(args, "minK", 4);
		runner.setKmerMin(minK);

		// Length of the largest K-mer in the K-mer models
		int maxK = Args.parseInteger(args, "maxK", 6);
		runner.setKmerMax(maxK);
		
		runner.setNumK();

		// K-mer models file / weights file
		String weights = Args.parseString(args, "weights", null);
		if (weights == null) {
			System.err.println("Provide weights file");
			System.exit(1);
		}
		runner.setKmerWeights(weights);
		
		// Now load K-mer models
		String motifFile = ap.getKeyValue("motiffile");
		runner.loadMotifsFromFile(motifFile);
		
		runner.execute();
	}

}
