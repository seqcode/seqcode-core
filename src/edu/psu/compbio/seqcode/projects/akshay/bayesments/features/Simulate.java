package edu.psu.compbio.seqcode.projects.akshay.bayesments.features;

import java.util.Random;

import edu.psu.compbio.seqcode.gse.utils.Pair;

public class Simulate {
	
	protected double[][] mean_matrix={{2.0,3.0,3.0,2.0,4.0},
			{2.0,3.0,2.0,4.0,1.0},
			{2.0,1.0,1.0,0.5,2.0}};
	protected double[][] sigma_matrix={{1.5,1.5,1.0,1.0,2.0},
			{1.5,1.5,0.5,2.0,0.5},
			{1.5,1.0,1.0,0.15,1.0}};
	protected int N=9000;
	protected int C=5;
	protected float[][] SimXc = new float[N][C-1];
	protected float[][] SimXf = new float[N][1];
	
	public Simulate() {
		// TODO Auto-generated constructor stub
	}
	
	public void simulate(){
		for(int c=0; c<C-1; c++){
			for(int i=0; i<3000; i++){ //Generating the FIRST out of the total 3 classes
				Pair<Float, Float> temp = this.doBoxMuller(mean_matrix[0][c], sigma_matrix[0][c]);
				SimXc[i][c] = temp.car();
				i++;
				SimXc[i][c] = temp.cdr();
			}	
		}
		for(int i=0; i<3000; i++){
			Pair<Float, Float> temp = this.doBoxMuller(mean_matrix[0][C-1], sigma_matrix[0][C-1]);
			SimXf[i][0] = temp.car();
			i++;
			SimXf[i][0] = temp.cdr();
		}
		
		for(int c=0; c<C-1; c++){
			for(int i=3000; i<5000; i++){ //Generating the SECOND out of the total 3 classes
				Pair<Float, Float> temp = this.doBoxMuller(mean_matrix[1][c], sigma_matrix[1][c]);
				SimXc[i][c] = temp.car();
				i++;
				SimXc[i][c] = temp.cdr();
			}	
		}
		for(int i=3000; i<5000; i++){
			Pair<Float, Float> temp = this.doBoxMuller(mean_matrix[1][C-1], sigma_matrix[1][C-1]);
			SimXf[i][0] = temp.car();
			i++;
			SimXf[i][0] = temp.cdr();
		}
		
		for(int c=0; c<C-1; c++){
			for(int i=5000; i<9000; i++){ //Generating the THIRD out of the total 3 classes
				Pair<Float, Float> temp = this.doBoxMuller(mean_matrix[2][c], sigma_matrix[2][c]);
				SimXc[i][c] = temp.car();
				i++;
				SimXc[i][c] = temp.cdr();
			}	
		}
		for(int i=5000; i<9000; i++){
			Pair<Float, Float> temp = this.doBoxMuller(mean_matrix[2][C-1], sigma_matrix[2][C-1]);
			SimXf[i][0] = temp.car();
			i++;
			SimXf[i][0] = temp.cdr();
		}
		
	}
	
	private Pair<Float,Float> doBoxMuller(double m, double s){
		float x1, x2, w, y1,y2;
		Random ran = new Random();
		do{
			x1 = (float) ( 2.0 * ran.nextFloat() - 1.0);
			x2 = (float) (2.0 * ran.nextDouble() - 1.0);
			w = x1 * x1 + x2 * x2;
		}while ( w >= 1.0 );
		w = (float) Math.sqrt( (-2.0 * Math.log(w) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;

		return new Pair<Float,Float>((float)(m+y1*s),(float)(m+y2*s));
	}
	
	//Accesors
	
	public float[][] getSimXc(){return this.SimXc;}
	public float[][] getSimXf(){return this.SimXf;}
	public int getNumTrainingEgs(){return this.N;}
	public int getNumChromCondition(){return this.C-1;}
	
	

}
