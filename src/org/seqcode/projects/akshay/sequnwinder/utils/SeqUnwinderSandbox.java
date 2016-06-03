package org.seqcode.projects.akshay.sequnwinder.utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.gse.utils.sequence.SequenceUtils;
import org.seqcode.projects.shaun.FreqMatrixImport;



public class SeqUnwinderSandbox {
	/** Learned k-mer weights */
	protected HashMap<String,double[]> weights = new HashMap<String,double[]>();
	protected List<String> kmerModelNames = new ArrayList<String>();
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	
	public int minK =4;
	public int maxK =5;
	public int numK = 0;
	
	//Settors
	public void setKmin(int mK){minK = mK;}
	public void setKmax(int mK){maxK = mK;}
	public void setNumK(){
		numK = 0;
		for(int k=minK; k<=maxK; k++ ){
			numK += (int)Math.pow(4, k);
		}
		
	}
	// Load freq matrices
	public void loadMotifsFromFile(String filename) throws NumberFormatException, IOException {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifs.addAll(motifImport.readTransfacMatricesAsFreqMatrices(filename));
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
            		weights.put(words[i], new double[numK]);
            	}
            }else{
            	if(!header){
            		System.err.println("Please provide a header in the K-mer weight file");
            		System.exit(1);
            	}
            	int ind = getKmerBaseInd(words[0]) + RegionFileUtilities.seq2int(words[0]);
            	for(int i = 1; i < words.length; i++ ){
            		weights.get(kmerModelNames.get(i-1))[ind] = Double.parseDouble(words[i]);
            	}
            }
           
        }reader.close();
	}
	
	public int getKmerBaseInd(String s){
		int baseInd = 0;
		for(int k=minK; k<s.length(); k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	
	public void scoreMotifs(){
		Set<String> modNames = weights.keySet();
		StringBuilder sb = new StringBuilder();
		sb.append("Motif");sb.append("\t");
		// Print header
		for(String s: modNames){
			sb.append(s);sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		for(WeightMatrix mot : motifs){
			sb.append(mot.getName());sb.append("\t");
			List<String> motKmers = WeightMatrix.getConsensusKmerList(mot, minK, maxK);
			
			// Print the k-mers
			System.out.println(mot.getName());
			for(String s : motKmers){
				System.out.println(s);
				System.out.println(SequenceUtils.reverseComplement(s));
			}
			
			for(String modName : modNames){
				double score=0.0;
				for(String s : motKmers){
					String revS = SequenceUtils.reverseComplement(s);
					int indS = RegionFileUtilities.seq2int(s);
					int indRevS = RegionFileUtilities.seq2int(revS);
					int KmerInd  = indS<indRevS ? indS : indRevS;
					int ind = getKmerBaseInd(s) + KmerInd;
					score = score + weights.get(modName)[ind];
				}
				sb.append(score);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}

		System.out.println(sb.toString());
	}
	
	
	public static void main(String[] args) throws NumberFormatException, IOException{
		ArgParser ap = new ArgParser(args);
		SeqUnwinderSandbox sbox = new SeqUnwinderSandbox();
		
		sbox.setKmin(Args.parseInteger(args, "minK", 4));
		sbox.setKmax(Args.parseInteger(args, "maxK", 5));
		sbox.setNumK();
		
		String weightsFile = ap.getKeyValue("weights");
		sbox.setKmerWeights(weightsFile);
		
		String motsFilename = ap.getKeyValue("motifs");
		sbox.loadMotifsFromFile(motsFilename);
		sbox.scoreMotifs();
		
	}
	
}
