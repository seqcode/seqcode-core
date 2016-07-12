package org.seqcode.genome.sequence.seqfunctions;

/**
 * Compares statistical properties of two CompositeSeqFunctions of same type
 * 
 * @author mahony
 *
 */
public class CompositeSeqFunctionComparison<F extends SeqFunction> {

	CompositeSeqFunction<F> signal;
	CompositeSeqFunction<F> background;
	SeqFunction function;
	int scoreWidth;
	double [][] pvalues;
	double [][] diffs;
	
	CompositeSeqFunctionComparison(CompositeSeqFunction<F> sig, CompositeSeqFunction<F> back){
		signal = sig;
		background = back;
		function = sig.getFunction();
		
		if(signal.getWidth()!=background.getWidth()){
			System.err.println("Trying to compare CompositeSeqFunctions of different widths!");
			System.exit(1);
		}else{
			scoreWidth=signal.getWidth();
			pvalues = new double[function.scoreDimension()][scoreWidth];
			diffs = new double[function.scoreDimension()][scoreWidth];
		}
		
	}
	
	public void execute(){
		//Diffs
		
		
		//Welch's T-tests
		
		//Multiple hypothesis correction
	}
}
