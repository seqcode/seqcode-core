package edu.psu.compbio.seqcode.projects.naomi.shapealign;

/**
 * SimilarityScore: Stores different similarity measures form the following paper
 * 
 * Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions: Sung-Hyuk Cha
 * 
 * @author naomi yamada
 *
 */

public class SimilarityScore {
	
	double x1,x2,y1,y2;
	double similarity_score;
	
	protected boolean euclidean_L2 = false;
	
	public SimilarityScore(double x_1, double x_2, double y_1, double y_2){
		
		this.x1 = x_1;
		this.x2 = x_2;
		this.y1 = y_1;
		this.y2 = y_2;		
		
	}
	
	// accessors
	public double getScore(){return similarity_score;}
	
	// setter
	public void setEuclideanL2(){euclidean_L2 = true;}
	
	public double computeScore(){
		
		double score = 0;
		if ( (x1 != x2) || (y1 != y2)){	
			if (euclidean_L2 == true){score = euclideanL2();}
			
		}

		return score;

	}
	
	
	public double euclideanL2(){
		
		double score = 1/(Math.pow(x1-x2, 2) + Math.pow(y1-y2, 2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;		
	}
	
	public double sorensen(){
		
		double score = (x1+x2+y1+y2)/(Math.abs(x1-x2) + Math.abs(y1-y2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;		
	}
	
	public double soergel(){
		
		double score = (Math.max(x1, x2) + Math.max(y1, y2))/(Math.abs(x1-x2) + Math.abs(y1-y2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;		
	}
	
	public double lorentzian(){
		
		double score = 1/(Math.log(1+Math.abs(x1-y2))+Math.log(1+Math.abs(y1-y2))) - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;		
	}
	
	public double PCE(){
		
		double pce =  (x1*x2 + y1*y2)/(Math.pow(x1, 2) + Math.pow(y1, 2) + Math.pow(x2, 2) + Math.pow(y2, 2) - ( x1*x2 + y1*y2 ));
				
		double score = pce - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;		
	}
	
	public double squared_chi(){

		double score =  1/(Math.pow(x1-x2, 2)/(x1+x2) + Math.pow(y1-y2, 2)/(y1+y2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;
	}
	
	public double divergence(){

		double denom =  2*(Math.pow(x1-x2, 2)/Math.pow(x1+x2,2) + Math.pow(y1-y2, 2)/Math.pow(y1+y2,2));
				
		double score = 1/denom - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;
	}
	
	public double clark(){

		double denom =  Math.sqrt(Math.pow(Math.abs(x1-x2)/(x1+x2), 2) + Math.pow(Math.abs(y1-y2)/(y1+y2), 2));
				
		double score = 1/denom - Math.abs(x1-x2)-Math.abs(y1-y2);
		
		return score;
	}
	

}
