package org.seqcode.genome.sequence.seqfunctions;

/**
 * Simply converts the sequence into a 4xL frequency matrix 
 * @author mahony
 *
 */
public class BaseFrequencyFunction implements SeqFunction{

	//Variables
	final int scoreDimension = 4;
	final int scoringOffset = 0;
	final int scoreWindowSize = 1;
	final boolean isBetweenNucs = false;
	final String[] labels = {"A", "C", "G", "T"};
	final String description = "Base frequencies";
	
	
	public double[][] score(String seq) throws SeqFunctionException {
		if(seq.length()<scoreWindowSize)
			throw new SeqFunctionException("Sequence too short for BaseFrequencyFunction");
		
		double [][] scores = new double[scoreDimension][seq.length()];
		String seqU = seq.toUpperCase();
		for(int i=0; i<seqU.length(); i++){
			char b = seqU.charAt(i);
			for(int x=0; x<scoreDimension; x++){ scores[x][i]=0;}
			
			switch(b){
				case 'A': scores[0][i]=1; break;
				case 'C': scores[1][i]=1; break;
				case 'G': scores[2][i]=1; break;
				case 'T': scores[3][i]=1; break;
				case 'N': default: scores[0][i]=0.25; scores[1][i]=0.25; scores[2][i]=0.25; scores[3][i]=0.25; break;
				case 'R' : scores[0][i]=0.5; scores[2][i]=0.5; break;
				case 'Y' : scores[3][i]=0.5; scores[1][i]=0.5; break;
				case 'M' : scores[0][i]=0.5; scores[1][i]=0.5; break;
				case 'K' : scores[2][i]=0.5; scores[3][i]=0.5; break;
				case 'S' : scores[1][i]=0.5; scores[2][i]=0.5; break;
				case 'W' : scores[0][i]=0.5; scores[3][i]=0.5; break;
			}
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
