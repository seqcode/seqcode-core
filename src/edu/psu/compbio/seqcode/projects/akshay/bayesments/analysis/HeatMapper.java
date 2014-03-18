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
	String Xtitle=null;
	String Ytitle=null;
	String[] xlab;
	Config conf;
	int nrow;
	int ncol;
	String chart_name;
	
	
	public HeatMapper(Config c, double[][] matrix, String Xtitle, String Ytitle, String[] xlab, String chart_name) {
		this.matrix = matrix;
		this.nrow = matrix.length;
		this.ncol = matrix[0].length;
		this.Xtitle = Xtitle;
		this.Ytitle = Ytitle;
		this.conf = c;
		this.xlab = xlab;
		this.chart_name = chart_name;
	}
	
	@SuppressWarnings("deprecation")
	public void plot(Color c, boolean log_scale){
		this.scaleMatrix();
		HeatChart map = new HeatChart(this.matrix);
		//map.setHighValueColour(new Color(221,20,20));
		map.setHighValueColour(c);
		//if(log_scale){
		//	map.setColourScale(map.SCALE_LOGARITHMIC);
		//}
		map.setChartMargin(600);
		map.setCellHeight(200);
		map.setCellWidth(200);
		
		map.setAxisLabelsFont(new Font("Ariel", Font.PLAIN, 55));
		
		if(Xtitle != null){
			map.setXAxisLabel(Xtitle);
		}
		
		if(xlab != null){
			map.setXValues(xlab);
		}
		map.setAxisValuesFont(new Font("Ariel", Font.PLAIN, 55));
		if(xlab != null){
			map.setShowXAxisValues(true);
		}else{
			map.setShowXAxisValues(false);
		}
		map.setShowYAxisValues(false);
		map.setBackgroundColour(new Color(0, 0, 0, 0));
		if(Ytitle != null){
			map.setYAxisLabel(Ytitle);
		}
		try {
			File f = new File(conf.getOutputImagesDir().getAbsoluteFile()+File.separator+this.chart_name+".png"); 
			map.saveToFile(f);
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
			int minind = BayesmentsSandbox.getMinindex(mu);
			for(int j=0; j<nrow; j++){
				matrix[j][c] = (mu[j]-mu[minind])/(mu[maxind]-mu[minind]);
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

	}
}
	
	

