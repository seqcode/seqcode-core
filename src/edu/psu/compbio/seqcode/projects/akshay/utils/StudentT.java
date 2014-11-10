package edu.psu.compbio.seqcode.projects.akshay.utils;

import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;

public class StudentT {
	public double[] vec1;
	public double[] vec2;
	public double pval;
	
	public StudentT(double[] vec1, double[] vec2) {
		this.vec1=vec1;
		this.vec2=vec2;
		
		double mean1 = StatUtil.mean(vec1);
		double mean2 = StatUtil.mean(vec2);
		
		double std1 = StatUtil.std(vec1);
		double std2 = StatUtil.std(vec2);
		
		double Tstat_num = Math.pow(mean1, 2.0) - Math.pow(mean2, 2.0);
		double Tstat_den = Math.sqrt((Math.pow(std1, 2.0)/vec1.length)+(Math.pow(std2, 2.0)/vec2.length));
		double Tstat = Tstat_num/Tstat_den;
		
		double v_num =  Math.pow((Math.pow(std1, 2.0)/vec1.length)+(Math.pow(std2, 2.0)/vec2.length),2);
		double v_den = (Math.pow((Math.pow(std1, 2.0)/vec1.length),2.0)/(vec1.length-1))+(Math.pow((Math.pow(std2, 2.0)/vec2.length),2.0)/(vec2.length-1));
		
		double v = v_num/v_den;
		pval = StatUtil.studentTPvalue(Tstat, v);
	}
	
	public double getPval(){return this.pval;}
	
	
	public static void main(String[] args){
		
		
	}
	
	
	

}
