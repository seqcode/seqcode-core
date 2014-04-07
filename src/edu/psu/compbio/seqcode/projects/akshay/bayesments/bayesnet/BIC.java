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
	public double[][] assignment;
	
	
	
	
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
		MAPassignment MLE = new MAPassignment(model, conf, null, this.numChromStates, this.numFacStates);
		MLE.execute(false);
		this.assignment = MLE.getMapAssignments();
		
	}
	
	public double calculateBicScore(){
		double bic = 0.0;
		bic = -2*this.getLikleHood()+((this.numChromStates*C+this.numFacStates*F+this.numChromStates*this.numFacStates)*Math.log(N));
		return bic;
	}
	
	public double calculateAicScore(){
		double aic  =  -2*this.getLikleHood()+(this.numChromStates*C+this.numFacStates*F+this.numChromStates*this.numFacStates);
		return aic;
	}
	
	private double getLikleHood(){
		double L=0;
		
		for(int i=0; i<N; i++){
			double P_x_c_b =0.0;
			double chromatinProduct = 1.0;
			double factorProduct = 1.0;
			double seqProduct = 1.0;
			for(int c=0; c<C; c++){
				NormalDistribution gaussian = new NormalDistribution(MUc[(int)assignment[i][0]][c],Math.pow(SIGMAc[(int)assignment[i][0]][c], 2.0));
				chromatinProduct = (c==0) ? gaussian.calcProbability((double) Xc[i][c]) : chromatinProduct*gaussian.calcProbability((double) Xc[i][c]);
			}
			for(int f=0; f<F; f++){
				NormalDistribution gaussian = new NormalDistribution(MUf[(int)assignment[i][1]][f],Math.pow(SIGMAf[(int)assignment[i][1]][f], 2.0));
				factorProduct = (f==0)? gaussian.calcProbability((double) Xf[i][f]) : factorProduct*gaussian.calcProbability((double) Xf[i][f]);
			}
			if(inSeqMode){
				for(int m=0; m<M; m++){
					NormalDistribution gaussian = new NormalDistribution(MUs[(int)assignment[i][0]][m],Math.pow(SIGMAs[(int)assignment[i][0]][m], 2.0));
					seqProduct = (m==0) ? gaussian.calcProbability((double) Xs[i][m]) : seqProduct*gaussian.calcProbability((double) Xs[i][m]);
				}
			}
			if(inSeqMode){
				P_x_c_b = Math.pow(chromatinProduct,conf.getChromWeight())*Math.pow(seqProduct,conf.getSeqWeight())*Bjk[(int)assignment[i][0]][(int)assignment[i][1]]*factorProduct;
			}else{
				P_x_c_b = Math.pow(chromatinProduct,conf.getChromWeight())*Bjk[(int)assignment[i][0]][(int)assignment[i][1]]*factorProduct;
			}
			L = L+Math.log(P_x_c_b);
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
