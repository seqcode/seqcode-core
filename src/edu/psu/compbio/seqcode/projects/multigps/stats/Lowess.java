package edu.psu.compbio.seqcode.projects.multigps.stats;

import Jama.Matrix;

/**
 * Lowess: This implementation of Lowess interpolation is taken from MIDAS (v2.22).
 * http://www.tm4.org/midas.html
 * 
 * References:
 * Quackenbush, J. Microarray data normalization and transformation. Nature Genetics. Vol.32 supplement pp496-501 (2002).
 * Yang, I.V. et al. Within the fold: assessing differential expression measures and reproducibility in microarray assays. Genome Biol. 3, research0062.1-0062.12 (2002).
 * Cleveland, W.S. Robust locally weighted regression and smoothing scatterplots. J. Amer. Stat. Assoc. 74, 829-836 (1979).
 *  
 * Copyright @ 2001-2003, The Institute for Genomic Research (TIGR).
 * All rights reserved.
 *
 * Lowess.java
 *
 * Created on December 10, 2001, 10:39 AM
 * @author  wliang
 * @version
 */
public class Lowess {
	private double bandwidthPct;
	private int windowSize;
	private static double[] xArray;
    private static double[] yArray;
    private static double[] yEst;
    private double[] xSub;
    private double[] ySub;
    private int subQueryPointIndex;
    private double[] weights;
    
    /**
     * Constructor: Note that xInput,yInput pairs should be sorted along the x-axis (increasing).
     * 
     * @param xInput
     * @param yInput
     * @param bandwidthPct : Reasonable values are 0.25 to 0.5
     */
    public Lowess(double[] xInput, double[] yInput, double bandwidthPct) {
        this.bandwidthPct = bandwidthPct;
        xArray = xInput;
        yArray = yInput;
        int totalDataCount = xArray.length;
        yEst = new double[totalDataCount];
        
        windowSize = new Double(bandwidthPct * totalDataCount).intValue();
        
        if (windowSize > 1){    //Issue Locfit2.0, 08-01-2002
            xSub = new double[windowSize];
            ySub = new double[windowSize];
            weights = new double[windowSize];
            
            //System.out.println("    ---- Lowess Bandwidth = " + windowSize);
            
            //Construct matrix
            Matrix X = new Matrix(windowSize, 2);
            Matrix WX = new Matrix(windowSize, 2);
            Matrix WY = new Matrix(windowSize, 1);
            
            for (int q = 0; q < totalDataCount; q++){
                
                getEstWindowSubArrayAndWeights(q, windowSize);
                
                //Construct matrix
                //Matrix X = new Matrix(windowSize, 2);
                //Matrix WX = new Matrix(windowSize, 2);
                //Matrix WY = new Matrix(windowSize, 1);
                
                for ( int row = 0; row < windowSize; row++ ) {
                    WX.set(row, 0, Math.sqrt(weights[row]));
                    WX.set(row, 1, Math.sqrt(weights[row]) * xSub[row]);
                    WY.set(row, 0, Math.sqrt(weights[row]) * ySub[row]);
                }
                
                //Apply first order estimation: local linear model
                Matrix WXT = WX.transpose();
                Matrix B = (WXT.times(WX)).inverse().times((WXT).times(WY));
                
                double B0 = (double)B.get(0, 0);
                double B1 = (double)B.get(1, 0);
                yEst[q] = B0 + B1 * xSub[subQueryPointIndex];
            }
        }else{ //windowSize = 0 or 1 //Issue Locfit2.0, 08-01-2002
            for (int q = 0; q < totalDataCount; q++){
                yEst[q] = 0;
            }
        }
    }
    
    /**
     * Estimate the y fit values at each x value
     * @param x
     * @return
     */
    public double[] estimateValues(double[] x){
    	double[] yVal = new double[x.length];
    	int totalDataCount = x.length;
    	if (windowSize > 1){ 
            xSub = new double[windowSize];
            ySub = new double[windowSize];
            weights = new double[windowSize];
            //Construct matrix
            Matrix WX = new Matrix(windowSize, 2);
            Matrix WY = new Matrix(windowSize, 1);
           
            for (int q = 0; q < totalDataCount; q++){
                getEstWindowSubArrayAndWeightsEV(x[q], windowSize);
                for ( int row = 0; row < windowSize; row++ ) {
                    WX.set(row, 0, Math.sqrt(weights[row]));
                    WX.set(row, 1, Math.sqrt(weights[row]) * xSub[row]);
                    WY.set(row, 0, Math.sqrt(weights[row]) * ySub[row]);
                }
                
                //Apply first order estimation: local linear model
                Matrix WXT = WX.transpose();
                Matrix B = (WXT.times(WX)).inverse().times((WXT).times(WY));
                
                double B0 = (double)B.get(0, 0);
                double B1 = (double)B.get(1, 0);
                yVal[q] = B0 + B1 * x[q];
            }
        }else{ //windowSize = 0 or 1 //Issue Locfit2.0, 08-01-2002
            for (int q = 0; q < totalDataCount; q++){
                yVal[q] = 0;
            }
        }
    	return yVal;
    }
    
    public double[] getXarray(){
        return xArray;
    }
    
    public double[] getYEst(){
        return yEst;
    }
    
    private void getEstWindowSubArrayAndWeights(int queryPointIndex, int window){
        
        int fullLength = xArray.length;
        int subStartIndex = 0;
        int windowSize = window;
        double maxDistance = 0;
        
        //Issue Locfit1.0, 08-01-2002
        //if ((queryPointIndex >= windowSize) && (queryPointIndex <= fullLength - windowSize)){//middle region data
        if ((queryPointIndex >= windowSize) && (queryPointIndex < fullLength - windowSize)){//middle region data
            for (int i = 0; i < windowSize; i++){
                if ( (xArray[queryPointIndex] - (xArray[queryPointIndex - windowSize + 1 + i])) <= (xArray[queryPointIndex + 1 + i] - xArray[queryPointIndex]) ){
                    subStartIndex = queryPointIndex - windowSize + 1 + i;
                    subQueryPointIndex = queryPointIndex - subStartIndex;
                    maxDistance = ( ( (xArray[queryPointIndex] - xArray[subStartIndex]) > (xArray[subStartIndex + windowSize - 1]) - xArray[queryPointIndex] ) ? (xArray[queryPointIndex] - xArray[subStartIndex]) : (xArray[subStartIndex + windowSize - 1] - xArray[queryPointIndex] ) );
                    break;
                }
            }
        }else if (queryPointIndex < windowSize){
            for (int i = 0; i < windowSize; i++){
                if ( (xArray[queryPointIndex] - xArray[i]) <= (xArray[windowSize + i] - xArray[queryPointIndex]) ){
                    subStartIndex = i;
                    subQueryPointIndex = queryPointIndex - subStartIndex;
                    maxDistance = ( (xArray[queryPointIndex] - xArray[subStartIndex]) > (xArray[subStartIndex + windowSize - 1] - xArray[queryPointIndex] ) ) ? (xArray[queryPointIndex] - xArray[subStartIndex]) : (xArray[subStartIndex + windowSize - 1] - xArray[queryPointIndex] );
                    break;
                }
            }
        }else{//right margin data
            for (int i = 0; i < windowSize; i++){
                if ( (xArray[fullLength - 1 - i] - xArray[queryPointIndex]) <= (xArray[queryPointIndex] - xArray[fullLength - windowSize - 1 - i]) ){
                    subStartIndex = fullLength - windowSize - i;
                    subQueryPointIndex = queryPointIndex - subStartIndex;
                    maxDistance = ( (xArray[queryPointIndex] - xArray[subStartIndex]) > (xArray[subStartIndex + windowSize - 1] - xArray[queryPointIndex] ) ) ? (xArray[queryPointIndex] - xArray[subStartIndex]) : (xArray[subStartIndex + windowSize - 1] - xArray[queryPointIndex] );
                    break;
                }
            }
        }
        
        System.arraycopy(xArray, subStartIndex, xSub, 0, windowSize);
        System.arraycopy(yArray, subStartIndex, ySub, 0, windowSize);
        
        //Calculate distances and find tricube weights
        for (int j = 0; j < windowSize; j ++) {
            double distance = Math.abs(xSub[subQueryPointIndex] - xSub[j]);
            weights[j] = (double)Math.pow(1.0 - Math.pow((distance / maxDistance), 3.0), 3.0);
        }
    }
    
    private void getEstWindowSubArrayAndWeightsEV(double queryX, int window){
        
        int fullLength = xArray.length;
        int subStartIndex = 0;
        int windowSize = window;
        double maxDistance = 0;
        
        int closestXPoint=0; double closestXDist=Double.MAX_VALUE; 
        for(int i=0; i<fullLength; i++){
        	if(Math.abs(xArray[i]-queryX)<closestXDist){
        		closestXDist = Math.abs(xArray[i]-queryX);
        		closestXPoint=i;
        	}
        }

        if ((closestXPoint >= windowSize) && (closestXPoint < fullLength - windowSize)){//middle region data
            for (int i = 0; i < windowSize; i++){
                if ( (xArray[closestXPoint] - (xArray[closestXPoint - windowSize + 1 + i])) <= (xArray[closestXPoint + 1 + i] - xArray[closestXPoint]) ){
                    subStartIndex = closestXPoint - windowSize + 1 + i;
                    subQueryPointIndex = closestXPoint - subStartIndex;
                    maxDistance = ( ( (queryX - xArray[subStartIndex]) > (xArray[subStartIndex + windowSize - 1]) - queryX ) ? (queryX - xArray[subStartIndex]) : (xArray[subStartIndex + windowSize - 1] - queryX ) );
                    break;
                }
            }
        }else if (closestXPoint < windowSize){
            for (int i = 0; i < windowSize; i++){
                if ( (xArray[closestXPoint] - xArray[i]) <= (xArray[windowSize + i] - xArray[closestXPoint]) ){
                    subStartIndex = i;
                    subQueryPointIndex = closestXPoint - subStartIndex;
                    maxDistance = ( (queryX - xArray[subStartIndex]) > (xArray[subStartIndex + windowSize - 1] - queryX ) ) ? (queryX - xArray[subStartIndex]) : (xArray[subStartIndex + windowSize - 1] - queryX );
                    break;
                }
            }
        }else{//right margin data
            for (int i = 0; i < windowSize; i++){
                if ( (xArray[fullLength - 1 - i] - xArray[closestXPoint]) <= (xArray[closestXPoint] - xArray[fullLength - windowSize - 1 - i]) ){
                    subStartIndex = fullLength - windowSize - i;
                    subQueryPointIndex = closestXPoint - subStartIndex;
                    maxDistance = ( (queryX - xArray[subStartIndex]) > (xArray[subStartIndex + windowSize - 1] - queryX ) ) ? (queryX - xArray[subStartIndex]) : (xArray[subStartIndex + windowSize - 1] - queryX );
                    break;
                }
            }
        }
        
        System.arraycopy(xArray, subStartIndex, xSub, 0, windowSize);
        System.arraycopy(yArray, subStartIndex, ySub, 0, windowSize);
        
        //Calculate distances and find tricube weights
        for (int j = 0; j < windowSize; j ++) {
            double distance = Math.abs(queryX - xSub[j]);
            weights[j] = (double)Math.pow(1.0 - Math.pow((distance / maxDistance), 3.0), 3.0);
        }
    }
}

