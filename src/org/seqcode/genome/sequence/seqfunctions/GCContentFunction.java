package org.seqcode.genome.sequence.seqfunctions;

/**
 * Scores the GC-content of a sequence using a sliding window of length W 
 * @author mahony
 *
 */
public class GCContentFunction implements SeqFunction{

	//Variables
	final int scoreDimension = 1;
	int scoringOffset = 0;
	int scoreWindowSize = 5;
	boolean isBetweenNucs = false;
	final String[] labels = {"GC"};
	final String description = "GC Content";
	
	public GCContentFunction(int W){
		scoreWindowSize = W;
		if(scoreWindowSize %2 ==0){
			scoringOffset = (scoreWindowSize/2)-1;
			isBetweenNucs=true;
		}else{
			scoringOffset = scoreWindowSize/2;
			isBetweenNucs=false;
		}
	}
	
	public double[][] score(String seq) throws SeqFunctionException {
		if(seq.length()<scoreWindowSize)
			throw new SeqFunctionException("Sequence too short for GCContentFunction");
		
		double [][] scores = new double[scoreDimension][seq.length()];
		double [] composition = new double[seq.length()]; 
		String seqU = seq.toUpperCase();
		for(int i=0; i<seqU.length(); i++){
			char b = seqU.charAt(i);
			scores[0][i]=0;
			
			switch(b){
				case 'A': composition[i]=0; break;
				case 'C': composition[i]=1; break;
				case 'G': composition[i]=1; break;
				case 'T': composition[i]=0; break;
				case 'N': default: composition[i]=0.5; break;
				case 'R' : composition[i]=0.5; break;
				case 'Y' : composition[i]=0.5; break;
				case 'M' : composition[i]=0.5; break;
				case 'K' : composition[i]=0.5; break;
				case 'S' : composition[i]=1; break;
				case 'W' : composition[i]=0; break;
			}
		}
		
		double win=(double)scoreWindowSize;
		for(int i=0; i<composition.length-scoreWindowSize+1; i++){
			double score=0;
			for(int w=0; w<scoreWindowSize; w++){
				score+=composition[i+w];
			}
			score/=win;
			scores[0][i+scoringOffset]=score;
		}
		
		return scores;
	}

	public int scoreDimension() {
		return scoreDimension;
	}

	public int scoringOffset() {
		return scoringOffset;
	}

	public int scoreWindowSize() {
		return scoreWindowSize;
	}

	public boolean isBetweenNucleotides() {
		return isBetweenNucs;
	}

	public String[] dimensionLabels() {
		return labels;
	}

	public String scoreDescription() {
		return description;
	}

}
