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
		String[] xlabs = new String[] {"K27me3","K4me1","K4me3","K27ac","DNase","K4me2","K9me3","K20me3","CTCF","OCT4","STAT3","SMAD1","SOX2","Cdx2-motif"};
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
		//double[][] matrix = new double[][] {{29.032, 13.7759, 10.089, 31.525, 265.071, 14.272},
		//		{4.401, 5.457, 3.183, 9.494, 33.9, 2.486},
		//		{3.152, 19.479, 2.592, 83.789, 309.339, 20.807},
		//		{1.801, 0.821, 1.808, 1.3977, 13.2268, 0.446}};
		//double [][] matrix = new double[][] {{0,0.053,0.015},
		//		{0.159,0.342,0.340},
		//		{0,0.442,0.323},
		//		{0.84,0.162,0.32}};
		double [][] matri = new double[][]{{1.724817459837947,1.4690214802134045,1.607873668625371,1.9992693728012991,3.834245952018254,0.9276876441588076,1.607873668625371,1.7396637595647115,0.9792997771557063,1.8542132603884676,1.8973831643951038,1.0103847991703656,1.8702755364669694,2.5083404166667527},
		{0.7201704303227637,0.5275041332147319,0.7024111041873385,0.5851642073790861,2.683529458002431,0.36013976054141794,0.7024111041873385,1.0176092371816727,0.7115209763408666,0.9249903388531917,1.1244803170649078,0.7090814465650671,1.0007126802966486,2.555150677717619},
		{1.7319984228567724,3.477079463315468,1.5198469829430195,4.393882748458835,5.978667044855343,3.3934128647064625,1.5198469829430195,1.718590053962126,1.4000013831249059,3.1412649166338316,2.8022104864934114,1.8924928758384263,3.335173450674261,2.0419008334291893}};
		
		double[][] matrix = new double[][] {{	0.9812481997990595,0.2747209353628731,0.25374008908058204},
				{0.018751178312337323,0.12384066063219165,0.7458420770572993},
				{0,0.6014384040049419,0}
		};
		
		
		HeatMapper map = new HeatMapper(matrix, "Experimental Track", "Chromatin State");
		map.plot();
	}
}
	
	

