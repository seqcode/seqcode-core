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
		map.setHighValueColour(new Color(20,20,221));
		map.setChartMargin(600);
		map.setCellHeight(200);
		map.setCellWidth(200);
		map.setAxisLabelsFont(new Font("Ariel", Font.PLAIN, 55));
	//	map.setXAxisLabel(Xtitle);
		//String[] xlabs = new String[] {"K27me3","K4me1","K4me3","K27ac","DNase","K4me2"};
		//map.setXValues(xlabs);
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
		//double[][] matrix = new double[][] {{29.032, 13.7759, 10.089, 31.525, 265.071, 14.272},
		//		{4.401, 5.457, 3.183, 9.494, 33.9, 2.486},
		//		{3.152, 19.479, 2.592, 83.789, 309.339, 20.807},
		//		{1.801, 0.821, 1.808, 1.3977, 13.2268, 0.446}};
		//double [][] matrix = new double[][] {{0,0.053,0.015},
		//		{0.159,0.342,0.340},
		//		{0,0.442,0.323},
		//		{0.84,0.162,0.32}};
		double [][] matri = new double[][]{{4.40118341884727,5.454400623480247,3.1828502751441996,9.489368630182394,33.88547988658285,2.4840755361408022},
		{3.1522251316091605,19.476689587319637,2.59224166155906,83.76512631956068,309.250191749674,20.80184934658802},
		{29.030666347292104,13.774674603349387,10.089122305375888,31.518710449288218,265.02065335492756,14.270024995137218},
		{1.8010107787097895,0.8209476253485581,1.8078185550881471,1.3971846878691898,13.224849491562479,0.4467358048566439}};
		
		double[][] matrix = new double[][] {{0.1599219574306244,0.34003915182638955,0.3424221327007531},
				{0,0.3241051228047087,0.44217886803934875},
				{0,0.015677709417432757,0.05303022441889957},
				{0.8400714404264896,0.3201780159514997,0.16236877484097206}
		};
		
		HeatMapper map = new HeatMapper(matrix, "Experimental Track", "Chromatin State");
		map.plot();
	}
}
	
	

