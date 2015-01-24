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
	
	public HeatMapper(double[][] matrix, String Xtitle, String Ytitle, String[] xlab, String chart_name) {
		this.matrix = matrix;
		this.nrow = matrix.length;
		this.ncol = matrix[0].length;
		this.Xtitle = Xtitle;
		this.Ytitle = Ytitle;
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
	
	public void plotlocal(Color c, boolean log_scale){
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
			File f = new File("/Users/akshaykakumanu/PSU/Bayesments/mu.png"); 
			System.out.println("akshay");
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
			if(mu[minind] < 0){
				for(int j=0; j<nrow; j++){
					matrix[j][c] = (mu[j]-mu[minind])/(mu[maxind]-mu[minind]);
				}
			}else{
				for(int j=0; j<nrow; j++){
					matrix[j][c] = mu[j]/mu[maxind];
				}
			}
		}
	}
	
	public static void main(String[] args){
		double[][] matrix = new double[][] {{1.4854176879668337,0.5460801406131417,0.49470346027397644,     0.5330971085147722 ,     0.5069919833339893},
				{2.979136730439021      , 2.8190744185409113  ,    1.5674775933583045  ,    2.4727405682637595  ,    0.8885798929603806},
				{1.6333962249543026     , 0.6580854562012127    ,  0.5056763562404272     , 1.8292784674348292    ,  0.5242121761834195},
				{0.29370014844108927   ,  0.4163582424009643 ,     0.45425452567466246  ,   0.7206624502795733   ,   1.5335507581505516},
				{2.469208864144808    ,   1.8080688054570715  ,    0.9402181609776454   ,   1.8258748161160494   ,   2.7023166788243413},
				{3.166922232056732      , 3.6866094866442642  ,    2.494527540013427    ,   2.6372088600535553    ,  2.490177717668567},
				{2.359953884071535   ,    4.929212820644236       ,5.487639746261503    ,   3.8913367722308734   ,   1.9452476395909675},
				{3.08228381719764     ,   4.371105630883926    ,   3.419140826438703   ,    3.5157309218612323    ,  1.1598669092389688},
				{1.0685877147721223  ,    0.49336997067365246  ,   0.4654807514942835   ,   1.725764642483823    ,   2.014225645438821},
				{1.869837960014813    ,   0.695711337614623   ,    0.5212556493857313   ,   0.523527977445202    ,   1.9642729514185078},
				{2.356152230741153     ,  1.2023264612394708   ,   0.617861667047022   ,    1.9980827085060142    ,  1.733535713680673}};
		//double [][] matrix = new double[][] {{0,0.053,0.015},
		//		{0.159,0.342,0.340},
		//		{0,0.442,0.323},
		//		{0.84,0.162,0.32}};
		double [][] matri = new double[][]{{0.3845485109629124,0.9347891676409756,0.022924703853385162},
		{0.21675228210416866,0.052686539795323166,0.9770752653135614},
		{0.3986992069329324,0.012524292563705321,0}};
		String[] xlabs = {"K4me1","K4me2","K4me3","K27ac","K27me3"};
		HeatMapper map = new HeatMapper(matrix, "marks", "states",xlabs, "mu");
		System.out.println("akshay");
		map.plotlocal(new Color(221,20,20),true);
		//double[][] matrix = new double[][] {{	0.9812481997990595,0.2747209353628731,0.25374008908058204},
		//		{0.018751178312337323,0.12384066063219165,0.7458420770572993},
		//		{0,0.6014384040049419,0}
		//};

	}
}
	
	

