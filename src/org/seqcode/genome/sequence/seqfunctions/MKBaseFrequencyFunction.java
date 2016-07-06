package org.seqcode.genome.sequence.seqfunctions;

/**
 * Converts the sequence into a 2xL matrix of aMino/Keto base frequencies 
 * @author mahony
 *
 */
public class MKBaseFrequencyFunction implements SeqFunction{

	//Variables
	final int scoreDimension = 2;
	final int scoringOffset = 0;
	final int scoreWindowSize = 1;
	final boolean isBetweenNucs = false;
	final String[] labels = {"M", "K"};
	final String description = "aMino/Keto base frequencies";
	
	
	public double[][] score(String seq) throws SeqFunctionException {
		if(seq.length()<scoreWindowSize)
			throw new SeqFunctionException("Sequence too short for MKBaseFrequencyFunction");
		
		double [][] scores = new double[scoreDimension][seq.length()];
		String seqU = seq.toUpperCase();
		for(int i=0; i<seqU.length(); i++){
			char b = seqU.charAt(i);
			for(int x=0; x<scoreDimension; x++){ scores[x][i]=0;}
			
			switch(b){
				case 'A': scores[0][i]=1; break;
				case 'C': scores[0][i]=1; break;
				case 'G': scores[1][i]=1; break;
				case 'T': scores[1][i]=1; break;
				case 'N': case 'S': case 'R': case 'W': case 'Y': default: scores[0][i]=0.5; scores[1][i]=0.5; break;
				case 'M' : scores[0][i]=1; break;
				case 'K' : scores[1][i]=1; break;
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
