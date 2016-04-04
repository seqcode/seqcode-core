package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.Stack;

/**
 * SmithWatermanAlignmentMatrix: build M and I matrix using SmithWaterman algorithms
 * 
 * @author naomi yamada
 *
 */

public class SmithWatermanAlignmentMatrix {
	
	protected float [][] regACounts;
	protected float [][] regBCounts;
	protected int window;
	
	protected Stack<Integer> traceBackTable = new Stack<Integer>();
	protected float maxScore;
	protected int alignStartXCoord;
	protected int alignStartYCoord;
	protected int alignEndXCoord;
	protected int alignEndYCoord;

	//constants for Smith-Waterman Algorithms
	final static float GAP_OPEN = 2;
	final static float GAP_EXT = 1;
	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;
	
	final static float MINIMUM_VALUE = Float.MIN_VALUE;
	
	public SmithWatermanAlignmentMatrix (float[][] regAarray, float[][] regBarray){
		
		regACounts = regAarray;
		regBCounts = regBarray;
		window = regAarray.length;
		
	}
	
	// accessors
	public float getMaxScore(){return maxScore;}
	public Stack<Integer> getTraceBack(){return traceBackTable;}
	public int getStartX(){return alignStartXCoord;}
	public int getStartY(){return alignStartYCoord;}
	public int getEndX(){return alignEndXCoord;}
	public int getEndY(){return alignEndYCoord;}
	
	//setters
	public void setMaxScore(float max){maxScore = max;}
	public void setTraceBack(Stack<Integer> traceBack){traceBackTable = traceBack;}
	public void setStartXCoord(int startXCoord){alignStartXCoord = startXCoord;}
	public void setStartYCoord(int startYCoord){alignStartYCoord = startYCoord;}
	public void setEndXCoord(int endXCoord){alignEndXCoord = endXCoord;}
	public void setEndYCoord(int endYCoord){alignEndYCoord = endYCoord;}
	
	public float computeScore(float aVal, float bVal){
		float score = (aVal + bVal)/2 - Math.abs(aVal - bVal);
		return score;		
	}
	
	public void buildMatrix(){
	
		//initialization of M[][] & I[][] & H[][] matrix
		float [][] M = new float [window+1][window+1];
		float [][] I = new float [window+1][window+1];
		for (int i = 0 ; i < window; i++){
			for (int j = 0 ; j < window; j++){
				M[i][j] = 0;
				I[i][j] = 0;
			}
		}
	
		//fill in M[i][j] & I[i][j] matrix
		for (int i = 1 ; i <window; i++){
			for (int j = 1 ; j <window ; j++){
				float mScore = computeScore(regACounts[i-1][0], regBCounts[i-1][0])
						+ computeScore(regACounts[i-1][1], regBCounts[i-1][1]);
			
				M[i][j] = Math.max(M[i-1][j-1] + mScore, I[i-1][j-1] + mScore);
			
				float temp_I[] = new float[3];
				temp_I[0] = M[i][j-1] - GAP_OPEN;
				temp_I[1] = I[i][j-1] - GAP_EXT;
				temp_I[2] = M[i-1][j] - GAP_OPEN;
				temp_I[3] = I[i-1][j] - GAP_EXT;
			
				float max_I = MINIMUM_VALUE;
				for (int k = 0 ; i < 4 ; i++){
					if (temp_I[k] > max_I){ max_I = temp_I[i];}
				}
				
				I[i][j] = max_I;
			}
		}

		// find the highest value
		float maxScore = 0;
		int x_coord = 0;
		int y_coord = 0;
	
		for (int i = Math.floorDiv(window, 2); i < window; i++){
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
		for (int j = Math.floorDiv(window, 2); j < window; j++){
			if (M[j][window] > maxScore){
				maxScore = M[j][window];
				x_coord = j;
				y_coord = window;	
			}
			if (I[j][window] > maxScore){
				maxScore = I[j][window];
				x_coord = j;
				y_coord = window;
			}
		}
		
		// back track to reconstruct the path
		float currentScore = maxScore;
		int current_x = x_coord;
		int current_y = y_coord;
		Stack<Integer> traceBack = new Stack<Integer>();
		
		int i = x_coord;
		int j = y_coord;
		
		while ( i != 0 || j != 0){

			float mScore = computeScore(regACounts[i-1][0], regBCounts[i-1][0])
					+ computeScore(regACounts[i-1][1], regBCounts[i-1][1]);
			
			// diagonal case
			if ( (M[i-1][j-1] + mScore == M[i][j]) || (I[i-1][j-1] + mScore == currentScore)){
				traceBack.push(DIAG);		
				currentScore = Math.max(M[i-1][j-1], I[i-1][j-1]);
			
			// left case
			}else if( (M[i][j-1]-GAP_OPEN == currentScore) || (I[i][j-1]-GAP_EXT == currentScore)){
				traceBack.push(LEFT);
				currentScore = Math.max(M[i][j-1]-1, I[i][j-1]);
			
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
