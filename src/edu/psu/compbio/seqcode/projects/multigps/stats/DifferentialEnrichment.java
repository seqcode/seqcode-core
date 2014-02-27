package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;

import edu.psu.compbio.seqcode.projects.shaun.viz.ScatterPlotConfigured;

/**
 * DifferentialEnrichment: parent class for testers of pairwise differential enrichment in count data
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public abstract class DifferentialEnrichment {

	protected final int MV_IMAGE_FITLINE_POINTS=200;
	protected String fileIDname="";
	
	public DifferentialEnrichment(){
		
	}
	
	/**
	 * All DifferentialEnrichment implementations require an execute methods that will take a CountsDataset, 
	 * process it for differential enrichment, and return a modified CountsDataset.
	 *  
	 * @param data : CountsDataset 
	 */
	public abstract CountsDataset execute(CountsDataset data);
	
	/**
	 * Save a mean-variance plot for a given condition
	 * @param xy
	 * @param fitLine
	 * @param conditionName
	 * @param rasterImage
	 */
	public void saveMeanVarPlot(Matrix xy, double[] yfit, String conditionName, boolean rasterImage){
		//Filename
		String fileName = conditionName+".mean-var";
		if(rasterImage)
			fileName = fileName+".png";
		else
			fileName = fileName+".svg";
		
		//Sub-sample fit-line and add (assuming xy is sorted, and yfit is coordinate matched with xy)
		Matrix fitxy;
		if(yfit.length>(MV_IMAGE_FITLINE_POINTS*2)){
			fitxy= new Matrix(MV_IMAGE_FITLINE_POINTS+1,2);
			for(int i=0; i<MV_IMAGE_FITLINE_POINTS; i++){
				int index = yfit.length/MV_IMAGE_FITLINE_POINTS*i;
				fitxy.set(i, 0, xy.get(index,0));
				fitxy.set(i, 1, yfit[index]);
			}
			fitxy.set(MV_IMAGE_FITLINE_POINTS, 0, xy.get(yfit.length-1,0));
			fitxy.set(MV_IMAGE_FITLINE_POINTS, 1, yfit[yfit.length-1]);
		}else{
			fitxy= new Matrix(yfit.length,2);
			for(int index=0; index<yfit.length; index++){
				fitxy.set(index, 0, xy.get(index,0));
				fitxy.set(index, 1, yfit[index]); 
			}
		}
		
		//Generate image
		ScatterPlotConfigured plotter = new ScatterPlotConfigured("MA plot");
		plotter.saveMeanVarPlot(xy, fitxy, fileName, rasterImage);
	}
	
	/**
	 * Print mean-variance values to files
	 * @param xy : 2D mean variance matrix
	 * @param conditionName : String
	 */
	public void printMeanVarData(Matrix xy, String conditionName){
		//Filename
		String fileName = conditionName+".mean-var.txt";
				
		//Print to file
		try {
			FileWriter fw = new FileWriter(fileName);
			fw.write("Point\tMean\tVar\n");
			for(int u=0; u<xy.getRowDimension(); u++)
				fw.write(u+"\t"+xy.get(u,0)+"\t"+xy.get(u,1)+"\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Set a file name identifier for the scripts that are produced by the differential enrichment implementation
	 */
	public void setFileIDname(String fid){fileIDname = fid;}
}
