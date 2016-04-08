package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.Stack;

/**
 * SmithWatermanAlignmentMatrix: builds M and I matrix using SmithWaterman algorithms
 * 
 * @author naomi yamada
 *
 */

public class SmithWatermanAlignmentMatrix {
	
	protected double [][] regACounts;
	protected double [][] regBCounts;
	protected int window;
	
	protected Stack<Integer> traceBackTable = new Stack<Integer>();
	protected double maxScore;
	protected int alignStartXCoord;
	protected int alignStartYCoord;
	protected int alignEndXCoord;
	protected int alignEndYCoord;

	//constants for Smith-Waterman Algorithms
//	final static double GAP_OPEN = 1;
//	final static double GAP_EXT = 0.5;
	
	final static double GAP_OPEN = 2;
	final static double GAP_EXT = 1;
	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;
	
	final static double MINIMUM_VALUE = -10000;
	
	public SmithWatermanAlignmentMatrix (double[][] regAarray, double[][] regBarray){
		
		regACounts = regAarray;
		regBCounts = regBarray;
		window = regAarray.length;
		
	}
	
	// accessors
	public double getMaxScore(){return maxScore;}
	public Stack<Integer> getTraceBack(){return traceBackTable;}
	public int getStartX(){return alignStartXCoord;}
	public int getStartY(){return alignStartYCoord;}
	public int getEndX(){return alignEndXCoord;}
	public int getEndY(){return alignEndYCoord;}
	
	//setters
	public void setMaxScore(double max){maxScore = max;}
	public void setTraceBack(Stack<Integer> traceBack){traceBackTable = traceBack;}
	public void setStartXCoord(int startXCoord){alignStartXCoord = startXCoord;}
	public void setStartYCoord(int startYCoord){alignStartYCoord = startYCoord;}
	public void setEndXCoord(int endXCoord){alignEndXCoord = endXCoord;}
	public void setEndYCoord(int endYCoord){alignEndYCoord = endYCoord;}
	
//	public double computeScore(double aVal, double bVal){
//		double score = (aVal + bVal)/2 - Math.abs(aVal - bVal);
//		return score;		
//	}
	
	public double computeScore(double aVal, double bVal){
	
	double score = MINIMUM_VALUE;
	if (aVal != bVal ){	
		score = 1/Math.sqrt(Math.abs(Math.pow(aVal,2)-Math.pow(bVal, 2)))*Math.pow((aVal+bVal)/2, 2)-Math.abs(aVal-bVal);
	}else{
		score = 20*Math.pow((aVal+bVal)/2, 2)-Math.abs(aVal-bVal);
	}		
	return score;
	}	
	
	// this needs to go with big GAP_OPEN and GAP_EXT
//	public double computeScore(double aVal, double bVal){
		
//		double score = MINIMUM_VALUE;
//		if (aVal != bVal ){	
//			score = 1/Math.sqrt(Math.pow(aVal-bVal, 2))*Math.pow((aVal+bVal)/2, 2)-Math.abs(aVal-bVal);
//		}else{
//			score = 100*Math.pow((aVal+bVal)/2, 2)-Math.abs(aVal-bVal);
//		}		
//		return score;
//	}	
	
	
	public void buildMatrix(){
		
		//initialization of M matrix
		double [][] M = new double [window+1][window+1];
		for (int i = 0 ; i <= window; i++){
			for (int j = 0 ; j <= window; j++)
				M[i][j] = 0;
		}	
		
		for (int i = 1 ; i <= window; i++){
			for (int j = 1 ; j <= window ; j++){
				double mScore = computeScore(regACounts[i-1][0], regBCounts[j-1][0])
						+ computeScore(regACounts[i-1][1], regBCounts[j-1][1]);
			
				double temp_M[] = new double[3];
				temp_M[0] = M[i-1][j-1] + mScore;
				temp_M[1] = M[i][j-1] - GAP_OPEN;
				temp_M[2] = M[i-1][j] - GAP_OPEN;
			
				double max_I = temp_M[0];
				for (int k = 1 ; k < 3 ; k++){

					if (temp_M[k] > max_I){ max_I = temp_M[k];}
				}
				M[i][j] = max_I;
			}
		}		

//		System.out.println("printing M matrix");
//		for (int i = 0; i <=window ; i++){
//			for (int j = 0; j <=window; j++){
//				System.out.print(M[i][j]+" ");
//			}
//			System.out.println();
//		}			

		// find the highest value
		double maxScore = MINIMUM_VALUE;
		int x_coord = 0;
		int y_coord = 0;
	
		for (int i = (int) Math.floor(window/2); i <= window; i++){
			if (M[i][window] > maxScore){
				maxScore = M[i][window];
				x_coord = i;
				y_coord = window;	
			}
		}
		for (int j = (int) Math.floor(window/2); j <= window; j++){
			if (M[window][j] > maxScore){
				maxScore = M[window][j];
				x_coord = window;
				y_coord = j;	
			}
		}
		
		System.out.println("max score is "+maxScore+" x coord is "+x_coord+" y coord is "+y_coord);
		
		// back track to reconstruct the path
		double currentScore = maxScore;
		int current_x = x_coord;
		int current_y = y_coord;
		Stack<Integer> traceBack = new Stack<Integer>();
		
		int i = x_coord;
		int j = y_coord;
		
		while ( i != 0 && j != 0){

			double mScore = computeScore(regACounts[i-1][0], regBCounts[j-1][0])
					+ computeScore(regACounts[i-1][1], regBCounts[j-1][1]);
			
			// diagonal case
			if ( M[i-1][j-1] + mScore == currentScore ){
				traceBack.push(DIAG);		
				currentScore = M[i-1][j-1];
			
			// left case
			}else if( M[i][j-1]-GAP_OPEN == currentScore ){
				traceBack.push(LEFT);
				currentScore = M[i][j-1];
			
			// right case
			}else{
				System.out.println("current score is "+currentScore);
				System.out.println("M value is"+M[i-1][j]);
				traceBack.push(UP);
				currentScore = M[i-1][j];
			}
			
			current_x = i;
			current_y = j;
			
			i--;
			j--;
		}
		
		setMaxScore(maxScore);
		setTraceBack(traceBack);
		setStartXCoord(current_x);
		setStartYCoord(current_y);
		setEndXCoord(x_coord);
		setEndYCoord(y_coord);
	}


	// builds M and I matrix; this is currently not working
	public void buildGappedMatrix(){
		
		//initialization of M[][] & I[][] & H[][] matrix
		double [][] M = new double [window+1][window+1];
		double [][] I = new double [window+1][window+1];
		for (int i = 0 ; i <= window; i++){
			for (int j = 0 ; j <= window; j++){
				M[i][j] = 0;
				I[i][j] = 0;
			}
		}	


		//fill in M[i][j] & I[i][j] matrix
		for (int i = 1 ; i <= window; i++){
			for (int j = 1 ; j <= window ; j++){
				double mScore = computeScore(regACounts[i-1][0], regBCounts[i-1][0])
						+ computeScore(regACounts[j-1][1], regBCounts[j-1][1]);
			
				M[i][j] = Math.max(M[i-1][j-1] + mScore, I[i-1][j-1] + mScore);
			
				double temp_I[] = new double[4];
				temp_I[0] = M[i][j-1] - GAP_OPEN;
				temp_I[1] = I[i][j-1] - GAP_EXT;
				temp_I[2] = M[i-1][j] - GAP_OPEN;
				temp_I[3] = I[i-1][j] - GAP_EXT;
			
				double max_I = MINIMUM_VALUE;
				for (int k = 0 ; k < 4 ; k++){
					if (temp_I[k] > max_I){ max_I = temp_I[k];}
				}
				
				I[i][j] = max_I;
			}
		}

		
		System.out.println("printing M matrix");
		for (int i = 0; i <=window ; i++){
			for (int j = 0; j <=window; j++){
				System.out.print(M[i][j]+" ");
			}
			System.out.println();
		}
		
		

		// find the highest value
		double maxScore = MINIMUM_VALUE;
		int x_coord = 0;
		int y_coord = 0;
	
		for (int i = (int) Math.floor(window/2); i <= window; i++){
			if (M[i][window] > maxScore){
				maxScore = M[i][window];
				x_coord = i;
				y_coord = window;	
			}
			if (I[i][window] > maxScore){
				maxScore = I[i][window];
				x_coord = i;
				y_coord = window;
			}
		}
		for (int j = (int) Math.floor(window/2); j <= window; j++){
			if (M[window][j] > maxScore){
				maxScore = M[window][j];
				x_coord = j;
				y_coord = window;	
			}
			if (I[window][j] > maxScore){
				maxScore = I[window][j];
				x_coord = j;
				y_coord = window;
			}
		}
		
		System.out.println("max score is "+maxScore+" x coord is "+x_coord+" y coord is "+y_coord);
		
		// back track to reconstruct the path
		double currentScore = maxScore;
		int current_x = x_coord;
		int current_y = y_coord;
		Stack<Integer> traceBack = new Stack<Integer>();
		
		int i = x_coord;
		int j = y_coord;
		
		while ( i != 0 && j != 0){

			double mScore = computeScore(regACounts[i-1][0], regBCounts[j-1][0])
					+ computeScore(regACounts[i-1][1], regBCounts[j-1][1]);
			
			// diagonal case
			if ( (M[i-1][j-1] + mScore == currentScore) || (I[i-1][j-1] + mScore == currentScore)){
				traceBack.push(DIAG);		
				currentScore = Math.max(M[i-1][j-1], I[i-1][j-1]);
			
			// left case
			}else if( (M[i][j-1]-GAP_OPEN == currentScore) || (I[i][j-1]-GAP_EXT == currentScore)){
				traceBack.push(LEFT);
				currentScore = Math.max(M[i][j-1], I[i][j-1]);
			
			// right case
			}else{
				traceBack.push(UP);
				currentScore = Math.max(M[i-1][j], I[i-1][j]);
			}
			
			current_x = i;
			current_y = j;
			
			i--;
			j--;
		}
		
		setMaxScore(maxScore);
		setTraceBack(traceBack);
		setStartXCoord(current_x);
		setStartYCoord(current_y);
		setEndXCoord(x_coord);
		setEndYCoord(y_coord);
	}
}