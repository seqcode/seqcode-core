package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class BIC {
	protected EMtrain model;
	protected float[][] Xc;
	protected float[][] Xf;
	protected double[][] Xs;
	
	protected double[][] MUc;
	protected double[][] MUf;
	protected double[][] MUs;
	
	protected double[][] SIGMAc;
	protected double[][] SIGMAf;
	protected double[][] SIGMAs;
	
	protected double[][] Bjk;
	
	protected int numChromStates;
	protected int numFacStates;
	
	protected Config conf;
	protected int N;
	protected int C;
	protected int F;
	protected int M;
	
	protected boolean inSeqMode;
	
	
	
	
	public BIC(Config conf, EMtrain model, int nChromStates, int nFacStates) {
		this.conf = conf;
		this.model = model;
		
		this.Xc = model.getXc();
		this.Xf = model.getXf();
		
		this.MUc = model.getMUc();
		this.MUf = model.getMUf();
		
		this.SIGMAc = model.getSIGMAc();
		this.SIGMAf = model.getSIGMAf();
		
		this.Bjk = model.getBjk();
		
		this.N = model.getNumTrainingEgs();
		this.C = model.getnumChromConds();
		this.F = model.getnumFacConds();
		
		this.inSeqMode = model.getSeqStateStatus();
		if(inSeqMode){
			this.setSeqParams();
		}
		this.numChromStates = nChromStates;
		this.numFacStates = nFacStates;
		
	}
	
	public double calculateBicScore(){
		double bic = 0.0;
		bic = -2*this.getLikleHood()+((this.numChromStates*C+this.numFacStates*F+this.numChromStates*this.numFacStates)*Math.log(N));
		return bic;
	}
	
	
	private double getLikleHood(){
		double L=0;
		
		for(int i=0; i<N; i++){
			double P_x=0.0;
			for(int j=0; j<this.numChromStates; j++){
				for(int k=0; k<this.numFacStates; k++){
					double P_x_c_b =0.0;
					double chromatinProduct = 1.0;
					double factorProduct = 1.0;
					double seqProduct = 1.0;
					for(int c=0; c<C; c++){
						NormalDistribution gaussian = new NormalDistribution(MUc[j][c],Math.pow(SIGMAc[j][c], 2.0));
						chromatinProduct = (c==0) ? gaussian.calcProbability((double) Xc[i][c]) : chromatinProduct*gaussian.calcProbability((double) Xc[i][c]);
					}
					for(int f=0; f<F; f++){
						NormalDistribution gaussian = new NormalDistribution(MUf[k][f],Math.pow(SIGMAf[k][f], 2.0));
						factorProduct = (f==0)? gaussian.calcLogProbability((double) Xf[i][f]) : factorProduct*gaussian.calcLogProbability((double) Xf[i][f]);
					}
					if(inSeqMode){
						for(int m=0; m<M; m++){
							NormalDistribution gaussian = new NormalDistribution(MUs[j][m],Math.pow(SIGMAs[j][m], 2.0));
							seqProduct = (m==0) ? gaussian.calcLogProbability((double) Xs[i][m]) : seqProduct*gaussian.calcLogProbability((double) Xs[i][m]);
						}
					}
					if(inSeqMode){
						P_x_c_b = chromatinProduct*seqProduct*Bjk[j][k]*factorProduct;
					}else{
						P_x_c_b = chromatinProduct*Bjk[j][k]*factorProduct;
					}
					P_x = P_x + P_x_c_b;
				}
			}
			L= L+Math.log(P_x);
		}
		return L;
		
	}
	
	
	private void setSeqParams(){
		this.Xs = model.getXs();
		this.MUs = model.getMUs();
		this.SIGMAs = model.getSIGMAs();
		this.M = model.getNumMotifs();
	}
	
}
