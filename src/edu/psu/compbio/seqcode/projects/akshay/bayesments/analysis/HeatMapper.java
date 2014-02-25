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
		String[] xlabs = new String[] {"K27me3","K4me1","K4me3","K27ac","DNase","K4me2","Cdx2-motif"};
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
		double [][] matri = new double[][]{{1.76150442496049,3.3261073705042556,1.466524328403386,4.025949670277766,5.500600575249554,3.0811686321305105,2.6567394896410756},
		{1.2801149781895813,0.8700851302741629,1.2164063944556835,1.2223294029367835,3.275982478775608,0.5099873625112245,2.934061798546805},
		{1.4686162118952928,1.708255759119876,1.5842543638468938,2.305979675560323,3.9669171220506696,1.2911155718953977,-1.4551907122765675}};
		
		double[][] matrix = new double[][] {{0.6311828793691772,0,0.1799676246134642},
				{0.3380894903205945,0.9198079157020873,0.6277287588766046},
				{0.030727630310349013,0.08014581094232466,0.19230361650994154}
		};
		
		
		HeatMapper map = new HeatMapper(matrix, "Experimental Track", "Chromatin State");
		map.plot();
	}
}
	
	

