package org.seqcode.genome.sequence.seqfunctions;

/**
 * Class of functions that score a sequence according to some metric.
 *  
 * Assumes that the function provides per-base (or between base) quantitative scores, albeit allowing for some offset. 
 *  
 * @author mahony
 *
 */
public interface SeqFunction {

	//Score a sequence, returning an array of score values
	public double[][] score(String seq) throws SeqFunctionException;
	
	//score dimensionality (per base)
	public int scoreDimension();
	//offset from initial base position to first score. 
	public int scoringOffset();
	//Window size is the window of sequence that each score is based on 
	public int scoreWindowSize();
	//the score is defined between nucleotide positions (i.e. first score defined at offset+0.5bp)
	public boolean isBetweenNucleotides();
	
	//Score max & min
	public double getMaxScore();
	public double getMinScore();
	
	//Labels on each dimension
	public String[] dimensionLabels();
	//Description of the score
	public String scoreDescription();
}
