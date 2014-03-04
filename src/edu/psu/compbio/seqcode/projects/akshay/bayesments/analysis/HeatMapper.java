package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;

import org.tc33.jheatchart.HeatChart;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;

public class HeatMapper {
	double[][] matrix;
	String Xtitle;
	String Ytitle;
	Config conf;
	int nrow;
	int ncol;
	
	
	public HeatMapper(double[][] matrix, String Xtitle, String Ytitle) {
		this.matrix = matrix;
		this.nrow = matrix.length;
		this.ncol = matrix[0].length;
		this.Xtitle = Xtitle;
		this.Ytitle = Ytitle;
		//this.conf = conf;
	}
	
	@SuppressWarnings("deprecation")
	public void plot(){
		this.scaleMatrix();
		HeatChart map = new HeatChart(this.matrix);
		map.setHighValueColour(new Color(221,20,20));
		map.setChartMargin(600);
		map.setCellHeight(200);
		map.setCellWidth(200);
		map.setAxisLabelsFont(new Font("Ariel", Font.PLAIN, 55));
		//map.setXAxisLabel(Xtitle);
		//String[] xlabs = new String[] {"K27me3","K4me1","K4me3","K27ac","DNase","K4me2","K9me3","K20me3","CTCF","OCT4","STAT3","SMAD1","SOX2","Cdx2-motif"};
		String[] xlabs = new String[] {"Sim-1","Sim-2","Sim-3","Sim-4"};
		map.setXValues(xlabs);
		map.setAxisValuesFont(new Font("Ariel", Font.PLAIN, 55));
		map.setShowYAxisValues(false);
		map.setShowXAxisValues(false);
		map.setBackgroundColour(new Color(0, 0, 0, 0));
		//map.setYAxisLabel(Ytitle);
		try {
			map.saveToFile(new File("/Users/akshaykakumanu/map.png"));
			//map.saveToFile(new File(conf.getOutputImagesDir().getAbsolutePath()+File.separator+"mu_heat.png"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private void scaleMatrix(){
		for(int c=0; c< ncol; c++){
			double[] mu = new double[nrow];
			for(int j=0; j<nrow; j++){
				mu[j] = matrix[j][c];
			}
			int maxind = BayesmentsSandbox.getMaxindex(mu);
			for(int j=0; j<nrow; j++){
				matrix[j][c] = mu[j]/mu[maxind];
			}
		}
	}
	
	public static void main(String[] args){
		double[][] matrix = new double[][] {{1.9937036105857127,3.1063555187858283,2.987112354682223,1.9806339284162375},
				{1.9918087985550605,1.0040493744787942,1.0211901276577968,0.5739006424128116},
			{1.9932481829956723,2.9657248091860255,1.990747548129314,5.133715938957541},};
		//double [][] matrix = new double[][] {{0,0.053,0.015},
		//		{0.159,0.342,0.340},
		//		{0,0.442,0.323},
		//		{0.84,0.162,0.32}};
		double [][] matri = new double[][]{{0.3845485109629124,0.9347891676409756,0.022924703853385162},
		{0.21675228210416866,0.052686539795323166,0.9770752653135614},
		{0.3986992069329324,0.012524292563705321,0}};
		
		//double[][] matrix = new double[][] {{	0.9812481997990595,0.2747209353628731,0.25374008908058204},
		//		{0.018751178312337323,0.12384066063219165,0.7458420770572993},
		//		{0,0.6014384040049419,0}
		//};
		
		
		HeatMapper map = new HeatMapper(matri, "Experimental Track", "Chromatin State");
		map.plot();
	}
}
	
	

