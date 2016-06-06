package org.seqcode.viz.scatter;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Paint;
import java.awt.geom.Ellipse2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.ui.TextAnchor;

import Jama.Matrix;

public class ScatterData {
	protected ArrayList<DefaultXYDataset> datasets = new ArrayList<DefaultXYDataset>();
	protected ArrayList<AnnotatedData> data = new ArrayList<AnnotatedData>();
	protected ArrayList<XYLineAndShapeRenderer> renderers = new ArrayList<XYLineAndShapeRenderer>();
    protected int datasetIndex = -1;

    public ScatterData(){
    	 
    }
    
    //Accessors
    public DefaultXYDataset getDataset(int i){if(datasetIndex>=0 && i>=0 && i<=datasetIndex){return datasets.get(i);}else{return null;}}
    public ArrayList<DefaultXYDataset> getDatasets(){return datasets;}
    public XYLineAndShapeRenderer getRenderer(int i){if(datasetIndex>=0 && i>=0 && i<=datasetIndex){return renderers.get(i);}else{return null;}}
    public ArrayList<XYLineAndShapeRenderer> getRenderers(){return renderers;}
    public int getDatasetIndex(){return datasetIndex;}
    
    /**
     * Add a dataset
     * @param dset
     * @parma drawLines: draw lines between points
     * @param drawShapes: draw the dots
     */
    public void addDataset(AnnotatedData dat){
    	data.add(dat);
    	DefaultXYDataset dset = new DefaultXYDataset();
    	dset.addSeries(dat.name, dat.values.transpose().getArray());
    	datasets.add(dset);
    	XYLineAndShapeRenderer currRend = new XYLineAndShapeRenderer(false, true);
    	currRend.setAutoPopulateSeriesFillPaint(false);
    	renderers.add(currRend);
    	datasetIndex++;
    }
    
    /**
     * Remove the last dataset added
     */
    public void removeLastDataset(){
    	data.remove(datasetIndex);
    	datasets.remove(datasetIndex);
    	renderers.remove(datasetIndex);
    	datasetIndex--;
    }
    
    /**
     * Return the minimum value in all datasets for dimension dim
     * @param dim
     * @return
     */
    public double getMin(int dim){
    	double min=Double.MAX_VALUE;
    	for(AnnotatedData d : data)
    		for(int x=0; x<d.values.getRowDimension(); x++)
    			if(d.values.get(x,dim) < min)
    				min=d.values.get(x,dim);
    	return min;
    }

    /**
     * Return the maximum value in all datasets for dimension dim
     * @param dim
     * @return
     */
    public double getMax(int dim){
    	double max=-Double.MAX_VALUE;
    	for(AnnotatedData d : data)
    		for(int x=0; x<d.values.getRowDimension(); x++)
    			if(d.values.get(x,dim) > max)
    				max=d.values.get(x,dim);
    	return max;
    }

    /**
     * Load a random dataset using createRandomDataset
     */
    public void loadRandomDataset(String name){
    	addDataset(createRandomDataset(name, 1000));
    }
    
    /**
     * Edit the parameters in the renderer
     * @param i
     * @param dotSize
     * @param p
     */
    public void editRenderer(int i, int dotSize, Paint p, boolean drawLines, boolean drawShapes){
    	XYLineAndShapeRenderer rend = renderers.get(i);
    	rend.setSeriesPaint(0, p);
    	rend.setSeriesShape(0, new Ellipse2D.Double(0, 0, dotSize, dotSize));
    	rend.setSeriesStroke(0, new BasicStroke(2));
    	rend.setSeriesLinesVisible(0, drawLines);
    	rend.setSeriesShapesVisible(0, drawShapes);
    }
    public void editRenderer(int i, int dotSize, Paint p){editRenderer(i, dotSize, p, false, true);}
    
    /**
     * Get annotated datapoints for a dataset
     * @param ID
     * @return
     */
    public List<XYTextAnnotation> getAnnotations(int ID, double xshift,double yshift, TextAnchor anchor, int textSize, int style, Color c){
    	List<XYTextAnnotation> annot = new ArrayList<XYTextAnnotation>();
    	for(int i=0; i<data.get(ID).annotations.length; i++){
    		XYTextAnnotation a = new XYTextAnnotation(data.get(ID).annotations[i], data.get(ID).annotationCoords.get(i,0)+xshift, data.get(ID).annotationCoords.get(i,0)+yshift);
    		a.setFont(new Font("Tahoma", style, textSize));
    		a.setTextAnchor(anchor);
    		a.setPaint(c);
    		annot.add(a);
    	}
    	return annot;
    }
    
    /**
     * Load a dataset from a 3-column file
     * @param f
     * @param skipLines
     */
    public void loadFileDataset(File f, String name, int skipLines){
    	try{
    		ArrayList<String> annots = new ArrayList<String>();
    		ArrayList<Double> Xdata=new ArrayList<Double>();
    		ArrayList<Double> Ydata=new ArrayList<Double>();
    		ArrayList<Double> Xac=new ArrayList<Double>();
    		ArrayList<Double> Yac=new ArrayList<Double>();
    		BufferedReader reader = new BufferedReader(new FileReader(f));
    		String line;
    		int count =0;
    		while ((line = reader.readLine()) != null) {
    			if(count>=skipLines){
	    			line = line.trim();
	    			String[] words = line.split("\\s+");
	    			if(words.length>=3){
	    				annots.add(words[0]);
	    				Xdata.add(new Double(words[1]));
	    				Ydata.add(new Double(words[2]));
	    				if(words.length>=5){
	    					Xac.add(new Double(words[3]));
		    				Yac.add(new Double(words[4]));
	    				}else{
	    					Xac.add(new Double(words[1]));
		    				Yac.add(new Double(words[2]));
	    				}
	    			}
    			}
    			count++;
    		}
    		int numPoints = Xdata.size();
    		AnnotatedData dat = new AnnotatedData(name, numPoints);
        	for(int i=0; i<numPoints; i++){
    			dat.annotations[i] = annots.get(i);
    			dat.values.set(i,0,Xdata.get(i));
    			dat.values.set(i,1,Ydata.get(i));
    			dat.annotationCoords.set(i,0,Xac.get(i));
    			dat.annotationCoords.set(i,1,Yac.get(i));
    		}
    		addDataset(dat);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    /**
     * Creates a dataset from a two dimensional array of datapoints.
     * Assumes that twoDimData is a 2D array (Mx2)
     * 
     * @param name  the series name.
     * 
     */
    public void loadDataset(String name, Matrix twoDimData) {
    	int numPoints = twoDimData.getRowDimension();
    	AnnotatedData dat = new AnnotatedData(name, numPoints);
    	for(int i=0; i<numPoints; i++){
    		double x = twoDimData.get(i, 0);
    		double y = twoDimData.get(i, 1);
    		dat.values.set(i,0,x);
    		dat.values.set(i,1,y);
    		dat.annotationCoords.set(i,0,x);
    		dat.annotationCoords.set(i,1,y);
    		dat.annotations[i] = ""+i;
    	}
        addDataset(dat);
    }
    
    /**
     * Creates a random dataset.
     * 
     * @param name  the series name.
     * 
     * @return The dataset.
     */
    private AnnotatedData createRandomDataset(String name, int numPoints) {
    	AnnotatedData dat = new AnnotatedData(name, numPoints);
    	Random generator = new Random();
    	for(int i=0; i<numPoints; i++){
    		double x = generator.nextGaussian();
    		double y = generator.nextGaussian();
    		dat.values.set(i,0,x);
    		dat.values.set(i,1,y);
    		dat.annotationCoords.set(i,0,x);
    		dat.annotationCoords.set(i,1,y);
    		dat.annotations[i] = "Rand"+i;
    	}
        return dat;
    }
    
    /**
     * AnnotatedData: datapoints with names
     * @author Shaun Mahony
     * @version	%I%, %G%
     */
    private class AnnotatedData{
    	public String name;
    	public Matrix values;
    	public String [] annotations;
    	public Matrix annotationCoords;
    	public AnnotatedData(String name, int numPoints){
    		this.name=name;
    		this.values = new Matrix(numPoints, 2);
    		this.annotationCoords = new Matrix(numPoints, 2);
    		this.annotations = new String[numPoints];
    	}
    }
}
