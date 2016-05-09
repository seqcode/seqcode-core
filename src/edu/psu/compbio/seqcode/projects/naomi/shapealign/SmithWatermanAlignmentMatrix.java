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
	
	static final double GAP_OPEN = 200;	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;	
	static final double MINIMUM_VALUE = -10000;
	
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
	

	public void buildMatrix(){
		
		//initialization of M matrix
		double [][] M = new double [window+1][window+1];
		for (int i = 0 ; i <= window; i++){
			for (int j = 0 ; j <= window; j++)
				M[i][j] = 0;
		}	
		
		for (int i = 1 ; i <= window; i++){
			for (int j = 1 ; j <= window ; j++){
				
				SimilarityScore similarity_s = new SimilarityScore(regACounts[i-1][0],regBCounts[j-1][0],regACounts[i-1][1],regBCounts[j-1][1]);
				
				similarity_s.setEuclideanL2();
				
				double mScore = similarity_s.computeScore();
			
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

			SimilarityScore similarity_s = new SimilarityScore(regACounts[i-1][0],regBCounts[j-1][0],regACounts[i-1][1],regBCounts[j-1][1]);
			
			similarity_s.setEuclideanL2();
			
			double mScore = similarity_s.computeScore();
			
			// diagonal case
			if ( M[i-1][j-1] + mScore == currentScore ){
				traceBack.push(DIAG);		
				currentScore = M[i-1][j-1];
				i--;
				j--;
			
			// left case
			}else if( M[i][j-1]-GAP_OPEN == currentScore ){
				traceBack.push(LEFT);
				currentScore = M[i][j-1];
				j--;
			
			// right case
			}else{
				System.out.println("current score is "+currentScore);
				System.out.println("M value is"+M[i-1][j]);
				traceBack.push(UP);
				currentScore = M[i-1][j];
				i--;
			}
			
			current_x = i;
			current_y = j;
		}
		
		setMaxScore(maxScore);
		setTraceBack(traceBack);
		setStartXCoord(current_x);
		setStartYCoord(current_y);
		setEndXCoord(x_coord);
		setEndYCoord(y_coord);
	}
}