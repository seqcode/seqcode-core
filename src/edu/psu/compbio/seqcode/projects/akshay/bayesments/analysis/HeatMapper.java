package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

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
	
	
	public HeatMapper(double[][] matrix, String Xtitle, String Ytitle, Config conf) {
		this.matrix = matrix;
		this.nrow = matrix.length;
		this.ncol = matrix[0].length;
		this.Xtitle = Xtitle;
		this.Ytitle = Ytitle;
		this.conf = conf;
	}
	
	public void plot(){
		this.scaleMatrix();
		HeatChart map = new HeatChart(this.matrix);
	}
	
	
	private void scaleMatrix(){
		for(int c=0; c< ncol; c++){
			double[] mu = new double[nrow];
			for(int j=0; j<nrow; j++){
				mu[j] = matrix[j][c];
			}
			int maxind = BayesmentsSandbox.getMaxindex(mu);
			for(int j=0; j<nrow; j++){
				mu[j] = mu[j]/mu[maxind];
			}
		}
	}
	
	
}
