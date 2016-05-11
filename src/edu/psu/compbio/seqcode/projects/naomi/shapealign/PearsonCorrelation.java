package edu.psu.compbio.seqcode.projects.naomi.shapealign;

public class PearsonCorrelation {
	
	protected double [][] normAarray;
	protected double [][] normBarray;
	protected double [][] revNormBarray;
	protected int window;
	
	protected boolean revAlign = false; 
	
	public PearsonCorrelation(double arrayA[][], double arrayB[][]){
		
		normAarray = arrayA;
		normBarray = arrayB;
		window = arrayB.length;
		
		revNormBarray = new double [window][2];
		
		for (int i = 0; i <= window; i++){
			for (int s = 0; s < 2 ;s++){
				revNormBarray[window-i][1-s] = arrayB[i][s];
			}
		}	
	}
	
	// accessor 
	public boolean isReverse(){return revAlign;}
	
	public double computePearsonCorr(){
		
		double [] catAarray = new double [2*(window+1)]; 
		double [] catBarray = new double [2*(window+1)]; 
		double [] catRevBarray = new double [2*(window+1)]; 
		
		// concatenate forward and reverse strands to make 1 dimensional arrays
		for (int i = 0 ; i <= window; i++){			
			catAarray[i] = normAarray[i][0];
			catAarray[i+window+1] = normAarray[i][1];
			
			catBarray[i] = normBarray[i][0];
			catBarray[i+window+1] = normBarray[i][1];
			
			catRevBarray[i] = revNormBarray[i][0];
			catRevBarray[i+window+1] = revNormBarray[i][1];	
		}

		double sumA = 0, sumB = 0;
		double aveA = 0, aveB = 0;		
		for (int i = 0 ; i < catAarray.length; i++){
			sumA += catAarray[i];
			sumB += catBarray[i];
		}
		aveA = sumA/catAarray.length;
		aveB = sumB/catBarray.length;
				
		// calculate Pearson correlation from forward and reverse pair
		double covf = 0, covr = 0;
		double varA = 0, varBf = 0, varBr = 0;
		double corrF = 0, corrR = 0;
		for (int i = 0; i < catAarray.length; i++){
			double ai = catAarray[i] - aveA;
			double bi_f = catBarray[i] - aveB;
			double bi_r = catRevBarray[i] - aveB;
			
			covf += covf + ai*bi_f;
			covr += covr + ai*bi_r;
			varA += ai*ai;
			varBf += bi_f*bi_f;
			varBr += bi_r*bi_r;
		}
		corrF = covf/(Math.sqrt(varA)*Math.sqrt(varBf));
		corrR = covr/(Math.sqrt(varA)*Math.sqrt(varBr));
		
		if (corrF > corrR){
			return corrF;
		}else{
			revAlign = true;
			return corrR;
		}
	}
}
