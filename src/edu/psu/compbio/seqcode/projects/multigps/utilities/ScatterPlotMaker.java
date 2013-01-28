package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberTickUnit;

import edu.psu.compbio.seqcode.projects.shaun.viz.ScatterPlot;

import Jama.Matrix;

/**
 * ScatterPlotterConfigured: Extension of ScatterPlot that implements a number of purpose-specific configurations as methods. 
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ScatterPlotMaker extends ScatterPlot{

	public ScatterPlotMaker(String title) {
		super(title);
	}
	
	/**
	 * Make an MA plot from one or two 2-D datasets (one will be highlighted in a different color) and save image.
	 * This configuration is used by deepseq.stats.Normalization classes
	 *  
	 * @param datapoints - 2D dataset (colored grey)
	 * @param datapoints_highlight - 2D dataset (colored blue), can be null
	 * @param yLine Double - data coordinates of line to be drawn parallel to x axis (used to show scaling line) 
	 * @param outFilename - String
	 * @param rasterImage - boolean
	 */
	public void saveMAplot(Matrix datapoints, Matrix datapoints_highlight, Double yLine, String outFilename, boolean rasterImage){
		this.setWidth(800);
		this.setHeight(800);
		this.setXLogScale(false);
		this.setYLogScale(false);
		this.addDataset("other", datapoints, new Color(75,75,75,80), 3);
		if(datapoints_highlight!=null)
			this.addDataset("highlighted", datapoints_highlight, new Color(0,0,255,80), 3);
		this.setXAxisLabel("A");
		this.setYAxisLabel("M");
		this.setXRangeFromData();
		this.setYRange(-10.5,10.5);

		//Set the tick units according to the range
		double xUpper = daxis.getRange().getUpperBound();
		double xLower = daxis.getRange().getLowerBound();
    	//if(daxis instanceof org.jfree.chart.axis.NumberAxis)
    	//	((NumberAxis)daxis).setTickUnit(new NumberTickUnit(5));
    	double yUpper = raxis.getRange().getUpperBound();
		double yLower = raxis.getRange().getLowerBound();
    	if(raxis instanceof org.jfree.chart.axis.NumberAxis)
    		((NumberAxis)raxis).setTickUnit(new NumberTickUnit(3));
		//Draw a line along y = yLine
    	if(yLine!=null){
    		XYLineAnnotation lineAnnot = new XYLineAnnotation(xLower, yLine, xUpper, yLine, new BasicStroke(1), new Color(0,0,0));
    		this.plot.addAnnotation(lineAnnot);
    	}
    	
    	try {
			this.saveImage(new File(outFilename), width, height, rasterImage);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Make an XY scatter plot from one or two 2-D datasets (one will be highlighted in a different color) and save image.
	 * This configuration is used by deepseq.stats.Normalization classes
	 *  
	 * @param datapoints - 2D dataset (colored grey)
	 * @param datapoints_highlight - 2D dataset (colored blue), can be null
	 * @param yLine Double - data coordinates of line to be drawn parallel to x axis (used to show scaling line) 
	 * @param outFilename - String
	 * @param rasterImage - boolean
	 */
	public void saveXYplot(Matrix datapoints, Matrix datapoints_highlight, String xName, String yName, String outFilename, boolean rasterImage){
		this.setWidth(800);
		this.setHeight(800);
		this.setXLogScale(false);
		this.setYLogScale(false);
		this.addDataset("other", datapoints, new Color(75,75,75,80), 3);
		if(datapoints_highlight!=null)
			this.addDataset("highlighted", datapoints_highlight, new Color(0,0,255,80), 3);
		this.setXAxisLabel(xName);
		this.setYAxisLabel(yName);
		this.setXRangeFromData();
		this.setYRangeFromData();
		
		//Set the tick units according to the range
		double xUpper = daxis.getRange().getUpperBound();
		double xLower = daxis.getRange().getLowerBound();
    	//if(daxis instanceof org.jfree.chart.axis.NumberAxis)
    	//	((NumberAxis)daxis).setTickUnit(new NumberTickUnit(5));
    	double yUpper = raxis.getRange().getUpperBound();
		double yLower = raxis.getRange().getLowerBound();
    	//if(raxis instanceof org.jfree.chart.axis.NumberAxis)
    	//	((NumberAxis)raxis).setTickUnit(new NumberTickUnit(3));
		
		try {
			this.saveImage(new File(outFilename), width, height, rasterImage);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
