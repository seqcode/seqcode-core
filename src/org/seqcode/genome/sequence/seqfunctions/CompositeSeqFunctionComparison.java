package org.seqcode.genome.sequence.seqfunctions;

import org.apache.commons.math3.stat.inference.TTest;

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
	int scoreWidth, zeroOffset;
	double [][] pvalues;
	double [][] diffs;
	
	public CompositeSeqFunctionComparison(CompositeSeqFunction<F> sig, CompositeSeqFunction<F> back){
		signal = sig;
		background = back;
		function = sig.getFunction();
		
		if(signal.getWidth()!=background.getWidth()){
			System.err.println("Trying to compare CompositeSeqFunctions of different widths!");
			System.exit(1);
		}else{
			scoreWidth=signal.getWidth();
			zeroOffset = signal.getZeroOffset();
			pvalues = new double[function.scoreDimension()][scoreWidth];
			diffs = new double[function.scoreDimension()][scoreWidth];
			for(int i=0; i<function.scoreDimension(); i++){
				for(int j=0; j<scoreWidth; j++){
					pvalues[i][j]=1.0;
					diffs[i][j]=0.0;
				}
			}
		}
	}
	
	//Accessors
	public SeqFunction getFunction(){return function;}
	public double[][] getPValues(){return pvalues;}
	public double[][] getDiffs(){return diffs;}
	public int getWidth(){return scoreWidth;}
	public int getZeroOffset(){return zeroOffset;}
	
	/**
	 * Calculate differences & p-values. 
	 */
	public void execute(){
		double [][] sigMeans = signal.getMeans();
		double [][] backMeans = background.getMeans();
		double [][] sigVars = signal.getVariances();
		double [][] backVars = background.getVariances();
		int sigNumSeq = signal.getNumSeqs();
		int backNumSeq = background.getNumSeqs();
		MyTTest test = new MyTTest();
		
		//Diffs
		for(int i=0; i<function.scoreDimension(); i++){
			for(int j=0; j<scoreWidth; j++){
				diffs[i][j]=(sigMeans[i][j]-backMeans[i][j])/(function.getMaxScore()-function.getMinScore());
			}
		}
		
		//Welch's T-tests (Bonferroni corrected)
		for(int i=0; i<function.scoreDimension(); i++){
			for(int j=0; j<scoreWidth; j++){
				pvalues[i][j]=test.welch(sigMeans[i][j], backMeans[i][j], sigVars[i][j], backVars[i][j], (double)sigNumSeq, (double)backNumSeq);
				pvalues[i][j] = Math.min(1.0,  pvalues[i][j]*scoreWidth);
			}
		}
		
	}
	
	public void printDiffs(){
		String[] labels = function.dimensionLabels();
		for(int i=0; i<function.scoreDimension(); i++){
			System.out.print(labels[i]);
			for(int j=0; j<scoreWidth; j++)
				System.out.print("\t"+diffs[i][j]);
			System.out.println("");
		}
	}

	public void printPVals(){
		String[] labels = function.dimensionLabels();
		for(int i=0; i<function.scoreDimension(); i++){
			System.out.print(labels[i]);
			for(int j=0; j<scoreWidth; j++)
				System.out.print("\t"+pvalues[i][j]);
			System.out.println("");
		}
	}

	private class MyTTest extends TTest{
		public MyTTest(){super();}
		
		public double welch(double m1, double m2, double v1, double v2, double n1, double n2){
			if(m1==0 && m2==0)
				return(1.0);
			else	
				return this.tTest(m1, m2, v1, v2, n1, n2);
		}
	}
}
