package org.seqcode.projects.naomi.shapealign;

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
	
	protected boolean linear = false;
	protected boolean euclidean_L2 = true;
	protected boolean sorensen = false;
	protected boolean soergel = false;
	protected boolean lorentzian = false;
	protected boolean pce = false;
	protected boolean squared_chi = false;
	protected boolean divergence = false;
	protected boolean clark = false;
	protected boolean AKL = false; // Average Kullback-Leibler
	protected boolean ALLR = false; // Average Log-likelihood ratio
	
	public SimilarityScore(){}		
	
	// accessors
	public double getScore(){return similarity_score;}
	
	// setter
	public void setLinear(){linear = true;}
	public void setEuclideanL2(){euclidean_L2 = true;}
	public void setSorensen(){sorensen = true;}
	public void setSoergel(){soergel = true;}
	public void setLorentzian(){lorentzian = true;}
	public void setPCE(){pce = true;}
	public void setSquaredChi(){squared_chi = true;}
	public void setDivergence(){divergence = true;}
	public void setClark(){clark = true;}
	public void setAKL(){AKL = true;}
	public void setALLR(){ALLR = true;}
	
	public double computeScore(double x_1, double x_2, double y_1, double y_2){		
		this.x1 = x_1;
		this.x2 = x_2;
		this.y1 = y_1;
		this.y2 = y_2;
				
		double score = 0;
		if ( (x1 != x2) || (y1 != y2)){	
			if (linear == true){score = linear();}
			else if (sorensen == true){score = sorensen();}
			else if (soergel == true){score = soergel();}
			else if (lorentzian == true){score = lorentzian();}
			else if (pce == true){score = PCE();}
			else if (squared_chi == true){score = squared_chi();}
			else if (divergence == true){score = divergence();}
			else if (clark == true){score = clark();}
			else if (AKL == true){score = AKL();}
			else if (ALLR == true){score = ALLR();}
			else{ score = euclideanL2();}			
		}
		return score;
	}
	
	protected double linear(){
		double score = (x1 + x2)/2 + (y1 + y2)/2 - (Math.abs(x1 - x2) + Math.abs(y1 - y2));
		return score;		
	}
		
	protected double euclideanL2(){		
		double score = 1 - Math.sqrt(Math.pow(x1-x2, 2) + Math.pow(y1-y2, 2)) - Math.abs(x1-x2)-Math.abs(y1-y2);		
		return score;		
	}
	
	protected double sorensen(){		
		double score = 1 - (Math.abs(x1-x2) + Math.abs(y1-y2))/(x1+x2+y1+y2) - Math.abs(x1-x2)-Math.abs(y1-y2);		
		return score;		
	}
	
	protected double soergel(){	
		double score = 1 - (Math.abs(x1-x2) + Math.abs(y1-y2))/(Math.max(x1, x2) + Math.max(y1, y2)) - Math.abs(x1-x2)-Math.abs(y1-y2);		
		return score;		
	}
	
	protected double lorentzian(){		
		double score = 1 - (Math.log(1+Math.abs(x1-y2))+Math.log(1+Math.abs(y1-y2))) - Math.abs(x1-x2)-Math.abs(y1-y2);		
		return score;		
	}
	
	protected double PCE(){		
		double pce =  (x1*x2 + y1*y2)/(Math.pow(x1, 2) + Math.pow(y1, 2) + Math.pow(x2, 2) + Math.pow(y2, 2) - ( x1*x2 + y1*y2 ));				
		double score = pce - Math.abs(x1-x2)-Math.abs(y1-y2);		
		return score;		
	}
	
	protected double squared_chi(){		
		double score = 0;		
		if (x1 == x2 && x1 == 0){
			score = 1 - Math.pow(y1-y2, 2)/(y1+y2) - Math.abs(y1-y2);
		}else if (y1 == y2 && y2 == 0){
			score = 1 - Math.pow(x1-x2, 2)/(x1+x2) - Math.abs(x1-x2);
		}else{
			score =  1 - (Math.pow(x1-x2, 2)/(x1+x2) + Math.pow(y1-y2, 2)/(y1+y2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		}	
		return score;
	}
	
	protected double divergence(){		
		double score = 0;
		if (x1 == x2 && x1 == 0){
			score = 2 - 2*Math.pow(y1-y2, 2)/Math.pow(y1+y2,2) - Math.abs(y1-y2);
		}else if (y1 == y2 && y2 == 0){
			score = 2 - 2*Math.pow(x1-x2, 2)/Math.pow(x1+x2,2) - Math.abs(x1-x2);
		}else{
			score = 2 - 2*(Math.pow(x1-x2, 2)/Math.pow(x1+x2,2) + Math.pow(y1-y2, 2)/Math.pow(y1+y2,2))- Math.abs(x1-x2)-Math.abs(y1-y2);
		}						
		return score;
	}
	
	protected double clark(){		
		double score = 0;
		if (x1 == x2 && x1 == 0){
			score = 1 - Math.sqrt(Math.pow(Math.abs(y1-y2)/(y1+y2), 2)) - Math.abs(y1-y2);
		}else if (y1 == y2 && y2 == 0){
			score = 1 - Math.sqrt(Math.pow(Math.abs(x1-x2)/(x1+x2), 2)) - Math.abs(x1-x2);
		}else{
			score = 1 - Math.sqrt(Math.pow(Math.abs(x1-x2)/(x1+x2), 2) + Math.pow(Math.abs(y1-y2)/(y1+y2), 2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		}			
		return score;
	}	
	
	protected double AKL(){
		double akl = 2 - (x1*Math.log(x1/x2) + x2*Math.log(x2/x1))/2 + (y1*Math.log(y1/y2) + y2*Math.log(y2/y1))/2;
		double score = akl - Math.abs(x1-x2)-Math.abs(y1-y2);
		return score;
	}
	
	protected double ALLR(){
		double score = 0;
		return score;
	}
}
