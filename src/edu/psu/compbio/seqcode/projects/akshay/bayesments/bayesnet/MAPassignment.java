package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class MAPassignment {
	
	//all the paramenters of the model
	
	public double[] PIj;
	public double[][] MUc;
	public double[][] MUf;
	public double[][] SIGMAc;
	public double[][] SIGMAf;
	public double[][] Bjk;
	public GenomicLocations trainingdata;
	public EMtrain model;
	public Config conf;
	public double[][] MapAssignment;
	
	public MAPassignment(EMtrain model, Config conf) {
		this.model = model;
		this.conf = conf;
		PIj = model.getPIj();
		MUc = model.getMUc();
		MUf = model.getMUf();
		SIGMAc = model.getSIGMAc();
		SIGMAf = model.getSIGMAf();
		Bjk = model.getBjk();
		trainingdata = model.getTrainingData();
		MapAssignment = new double[trainingdata.getNumTrainingExamples()][2];
	}
	
	public void execute(){
		int N = trainingdata.getNumTrainingExamples();
		int numChromStates = conf.getNumChrmStates();
		int numFacState = conf.getNumFacStates();
		int C = trainingdata.getNumChromatinCons();
		int F = trainingdata.getNumFacCons();
		float[][] Xc = trainingdata.getChromatinCounts();
		float[][] Xf = trainingdata.getFactorCounts();
		
		for(int i=0; i<N; i++){
			double[] assignment = new double[2];
			double maxLiklehood = 0.0;
			for(int j=0; j<numChromStates; j++){
				for(int k=0; k<numFacState; k++){
					double liklehood = 0.0;
					double chromGausssianProd=0.0;
					double facGaussianProd = 0.0;
					for(int c=0; c<C; c++){
						NormalDistribution gaussian = new NormalDistribution(MUc[j][c],Math.pow(SIGMAc[j][c], 2.0));
						chromGausssianProd = (c==0 ? gaussian.calcProbability((double) Xc[i][c]): chromGausssianProd* gaussian.calcProbability((double) Xc[i][c]));
					}
					for(int f=0; f< F; f++){
						NormalDistribution gaussian = new NormalDistribution(MUf[k][f],Math.pow(SIGMAf[k][f], 2.0));
						facGaussianProd = (f == 0 ? gaussian.calcProbability((double) Xf[i][f]): facGaussianProd* gaussian.calcProbability((double) Xf[i][f]));
					}
					chromGausssianProd = ( Double.isNaN(chromGausssianProd)) ? 0.0 : chromGausssianProd;
					facGaussianProd = (Double.isNaN(facGaussianProd)) ? 0.0: facGaussianProd;
					liklehood = PIj[j]*chromGausssianProd*Bjk[j][k]*facGaussianProd;
					liklehood = Double.isNaN(liklehood) ? 0 : liklehood;
					if(liklehood > maxLiklehood){
						maxLiklehood = liklehood;
						assignment[0] = j;
						assignment[1] = k;
					}
				}
			}
		}
	}
	
	

}
