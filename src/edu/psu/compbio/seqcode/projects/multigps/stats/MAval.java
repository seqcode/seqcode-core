package edu.psu.compbio.seqcode.projects.multigps.stats;

/**
 * MAval: values of an M-A plot
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class MAval{
	public double x=0;
	public double y=0;
	public double log_r=0;
	public double log_g=0;
	public double M=0;
	public double A=0;
	public double w=0;
	protected static double LOG_2 = Math.log(2.0);
	
	public MAval(double i, double j, double total_i, double total_j){
		x=i; y=j;
		if(x>0 && y>0){
			log_r = Math.log(i/total_i)/LOG_2;
			log_g = Math.log(j/total_j)/LOG_2;
			M = log_r - log_g;
			A = 0.5* (log_r + log_g);
			w = ((total_i-i)/(total_i*i))+((total_j-j)/(total_j*j));
		}
	}
	public int compareByM(MAval k){
		if(this.M < k.M){return -1;}
		else if(this.M > k.M){return 1;}
		else{return 0;}
	}
	public int compareByA(MAval k){
		if(this.A < k.A){return -1;}
		else if(this.A > k.A){return 1;}
		else{return 0;}
	}
}
