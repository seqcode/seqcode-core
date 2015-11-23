package edu.psu.compbio.seqcode.projects.akshay.MultiSeq;

import java.util.ArrayList;

import weka.core.Optimization;
import weka.core.Utils;
import edu.psu.compbio.seqcode.projects.akshay.MultiSeq.SeqUnwinder.ClassRelationStructure;
import edu.psu.compbio.seqcode.projects.akshay.MultiSeq.SeqUnwinder.ClassRelationStructure.Node;

public class AdmmFramework {
	
	public final double ADMM_PHO = 1.2; 
	public final double ADMM_ALPHA = 1.5; // Over-relaxation parameter (helps faster convergence)
	public final double ADMM_ABSTOL = 1E-4; // Tolerance parameters (for the stopping criteria)
	public final double ADMM_RELTOL = 1E-2;
	
	public int maxItr = 1000; // maximum allowed iterations for the ADMM algorithm
	
	public int BGFS_maxIts=-1;
	
	public double[] sm_x; // predictor weights of all nodes
	public double[] x; // predictor weights of only the leaf nodes
	public double[] z; // for the Z-step 
	public double[] u; // Augmented Lagrange dual variables
	
	public double[][] data;
	
	/** Weights of the instances */
	protected double[] weights;
	
	/** Instance class membership */
	protected int[] cls;
	
	/** Flowing are a bunch of SeqUnwinder model parameters */
	protected int numPredictors;
	protected int numClasses;
	protected ClassRelationStructure classStructure;
	protected int numNodes;
	protected double ridge;
	
	// To clear memory
	public void clearADMM(){
		data=null;
		sm_x=null;
		x=null;
		z=null;
		u=null;
	}
	
	//Settors
	public void setBGFSmaxItrs(int m){BGFS_maxIts=m;}
	public void setADMMmaxItrs(int m){maxItr = m;}
	public void setInstanceWeights(double[] w){weights=w;}
	public void setClsMembership(int[] c){cls=c;}
	public void setNumPredictors(int p){numPredictors=p;}
	public void setNumClasses(int c){numClasses = c;}
	public void setClassStructure(ClassRelationStructure rel){classStructure = rel; setNumNodes(rel.allNodes.size());}
	public void setNumNodes(int n){numNodes = n;}
	public void setRidge(double r){ridge = r;}
	
	//gettors
	public double[] getX(){return x;}
	public double[] getZ(){return z;}
	public double[] getU(){return u;}
	public double[] getsmX(){return sm_x;}
	
	// Initialize
	public void initZandU(){
		int dim = numPredictors+1;
		z= new double[numNodes*numNodes*dim];
		u= new double[numNodes*numNodes*dim];
	}
	
	
	
	public AdmmFramework(double[] xinit, double[] sm_xinit, double[][] d) {
		x = xinit;
		sm_x = sm_xinit;
		data=d;
	}
	
	
	
	public void execute() throws Exception{ // The main ADMM algorithm
		
		int dim = numPredictors +1;
		
		for(int ad=0; ad < maxItr; ad++){
			// x-update
			updateX();
		
			//z-update
			double[] zold = z;
			double[] xRelaxed = new double[z.length];
			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() > 0){
					for(int pid: n.parents){
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						int pOffset = pid*dim;
						for(int w=0; w<dim; w++){
							xRelaxed[zOffset+w] = ADMM_ALPHA*(x[nOffset+w]-x[pOffset+w])+(1-ADMM_ALPHA)*z[zOffset+w]+u[zOffset+w];
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=0; w<dim; w++){
						xRelaxed[zOffset+w] = ADMM_ALPHA*(x[nOffset+w])+(1-ADMM_ALPHA)*z[zOffset+w]+u[zOffset+w];
					}
				}
			}
		
			doSkrinkage(xRelaxed,ADMM_PHO/(2*ridge));
		
			//u-update
			for(Node n : classStructure.leafs){
				if(n.parents.size() > 0){
					for(int pid: n.parents){
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						for(int w=0; w<dim; w++){
							u[zOffset+w] = xRelaxed[zOffset+w]-z[zOffset+w];
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=0; w<dim; w++){
						u[zOffset+w] = xRelaxed[zOffset+w]-z[zOffset+w];
					}
				}
			}
		
			//Check for stopping criteria
			double[] primal_tol = new double[numNodes*numNodes];
			double[] dual_tol = new double[numNodes*numNodes];
		
			double[] primals = new double[numNodes*numNodes] ;
			double[] duals = new double[numNodes*numNodes];
		
			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() > 0){
					for(int pid: n.parents){
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						int pOffset = pid*dim;
						double xnorm=0;
						double znorm=0;
						double unorm=0;
						for(int w=0; w<dim; w++){
							xnorm +=  (x[nOffset+w]-x[pOffset+w])*(x[nOffset+w]-x[pOffset+w]);
							znorm += z[zOffset+w]*z[zOffset+w];
							unorm += ADMM_PHO*u[zOffset+w]*ADMM_PHO*u[zOffset+w];
						}
						xnorm = Math.sqrt(xnorm);
						znorm = Math.sqrt(znorm);
						unorm = Math.sqrt(unorm);
						primal_tol[n.nodeIndex*numNodes+pid] = Math.sqrt(numPredictors)*ADMM_ABSTOL+ADMM_RELTOL*Math.max(xnorm, znorm);
						dual_tol[n.nodeIndex*numNodes+pid] = Math.sqrt(numPredictors)*ADMM_ABSTOL + ADMM_RELTOL*unorm;
					
						double deltaZnorm=0;
						double primalResidueNorm = 0;
						for(int w=0; w<dim; w++){
							deltaZnorm +=  ADMM_PHO*(z[zOffset+w] - zold[zOffset+w])*ADMM_PHO*(z[zOffset+w] - zold[zOffset+w]);
							primalResidueNorm += (x[nOffset+w]-x[pOffset+w]-z[zOffset+w])*(x[nOffset+w]-x[pOffset+w]-z[zOffset+w]);
						}
						deltaZnorm = Math.sqrt(deltaZnorm);
						primalResidueNorm = Math.sqrt(primalResidueNorm);
						duals[n.nodeIndex*numNodes+pid] = deltaZnorm;
						primals[n.nodeIndex*numNodes+pid] = primalResidueNorm;
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					double xnorm=0;
					double znorm=0;
					double unorm=0;
					for(int w=0; w<dim; w++){
						xnorm +=  (x[nOffset+w])*(x[nOffset+w]);
						znorm += z[zOffset+w]*z[zOffset+w];
						unorm += ADMM_PHO*u[zOffset+w]*ADMM_PHO*u[zOffset+w];
					}
					xnorm = Math.sqrt(xnorm);
					znorm = Math.sqrt(znorm);
					unorm = Math.sqrt(unorm);
					primal_tol[n.nodeIndex*numNodes+n.nodeIndex] = Math.sqrt(numPredictors)*ADMM_ABSTOL+ADMM_RELTOL*Math.max(xnorm, znorm);
					dual_tol[n.nodeIndex*numNodes+n.nodeIndex] = Math.sqrt(numPredictors)*ADMM_ABSTOL + ADMM_RELTOL*unorm;
					double deltaZnorm=0;
					double primalResidueNorm = 0;
					for(int w=0; w<dim; w++){
						deltaZnorm +=  ADMM_PHO*(z[zOffset+w] - zold[zOffset+w])*ADMM_PHO*(z[zOffset+w] - zold[zOffset+w]);
						primalResidueNorm += (x[nOffset+w]-z[zOffset+w])*(x[nOffset+w]-z[zOffset+w]);
					}
					deltaZnorm = Math.sqrt(deltaZnorm);
					primalResidueNorm = Math.sqrt(primalResidueNorm);
					duals[n.nodeIndex*numNodes+n.nodeIndex] = deltaZnorm;
					primals[n.nodeIndex*numNodes+n.nodeIndex] = primalResidueNorm;
				}
			}
		
			boolean converged = true;
			for(int d=0; d<duals.length; d++ ){
				if(duals[d] > dual_tol[d])
					converged=false;
				if(primals[d] > primal_tol[d])
					converged=false;
				if(!converged)
					break;
			}
		
			if(converged)
				break;
		}
	}

	
	public void updateX() throws Exception{
		OptObject oO = new OptObject();
		Optimization opt = new OptEng(oO);
		opt.setMaxIteration(BGFS_maxIts);
		
		double[][] b = new double[2][x.length]; // Boundary constraints, N/A here
		for (int p = 0; p < numClasses; p++) {
			int offset = p * (numClasses + 1);
			b[0][offset] = Double.NaN;
			b[1][offset] = Double.NaN;
			for (int q = 1; q <= numPredictors; q++) {
				b[0][offset + q] = Double.NaN;
				b[1][offset + q] = Double.NaN;
			}
		}
		if (BGFS_maxIts == -1) { // Search until convergence
			x = opt.findArgmin(x, b);
			while (x == null) {
				x = opt.getVarbValues();
				x = opt.findArgmin(x, b);
			}
		} else {
			opt.setMaxIteration(BGFS_maxIts);
			x = opt.findArgmin(x, b);
			if (x == null) {
				x = opt.getVarbValues();
			}
		}
		// Just to make sure all the internal nodes have been updated
		oO.updateInternalNodes(x);
	}

	
	/**
	 * Essentially the Z-update step in the ADMM framework
	 * @param w
	 * @param pho
	 * @return
	 */
	public void doSkrinkage(double[] xrel, double pho){ 
		for(int i=0; i<xrel.length; i++){
			z[i] = xrel[i] - Math.signum(xrel[i])*Math.min(pho, xrel[i]);
		}
	}
	
	public class OptObject {
		
		
		public OptObject() {
		}
		
		/** 
		 * Fix the leaf node values
		 * Updates all odd numbered layers
		 * Update all even numbered layers except the leaf 
		 * (The gradient for the leaf layer is computed in the evaluateGradient code)
		 * @param x // Leaf layer nodes
		 */
		protected void updateInternalNodes(double[] currx){
			int dim = numPredictors+1;
			
			// Copy the current x(leaf params) to x
			for(Node l: classStructure.leafs){
				int offset = l.nodeIndex*dim;
				for(int w=0; w<dim;w++){
					sm_x[offset+w] = currx[offset+w];
				}
			}
			
			//First update all odd-numbered layrers
			for(int l=1; l< classStructure.numLayers;l+=2){
				for(Node n : classStructure.layers.get(l)){// Get nodes in this layer
					int offset = n.nodeIndex*dim;
					double den = n.parents.size()+n.children.size();
					for(int w=0; w<dim; w++){
						double num = 0.0;
						for(int pind: n.parents){
							num = num+sm_x[pind*dim+w];
						}
						for(int cind: n.children){
							num=num+sm_x[cind*dim+w];
						}
						sm_x[offset+w] = num/den;
					}
				}
			}
			
			// Now update all the even numbered layers except the leaf layer
			for(int l=2; l< classStructure.numLayers; l+=2){
				for(Node n:  classStructure.layers.get(l)){
					int offset = n.nodeIndex*dim;
					double den = n.parents.size()+n.children.size();
					for(int w=0; w<dim; w++){
						double num = 0;
						for(int pind : n.parents){
							num=num+sm_x[pind*dim+w];
						}
						for(int cind : n.children){
							num=num+sm_x[cind*dim+w];
						}
						sm_x[offset+w] = num/den;
					}
				}
			}
		}
		
		public double[] evaluateGradient(double[] currx){
			// Update the internal nodes first
			updateInternalNodes(currx);
			double[] grad = new double[currx.length];
			int dim = numPredictors + 1; // Number of variables per class

			for (int i = 0; i < cls.length; i++) { // ith instance
				double[] num = new double[numClasses]; // numerator of
		                                                     // [-log(1+sum(exp))]'
		        int index;
		        for (int offset = 0; offset < numClasses; offset++) { // Which
		                                                                    // part of 
		        	double exp = 0.0;
		        	index = offset * dim;
		        	for (int j = 0; j < dim; j++) {
		        		exp += data[i][j]*sm_x[index + j];
		        	}
		        	num[offset] = exp;
		        }

		        double max = num[Utils.maxIndex(num)];
		       
		        double denom=0.0;
		        for (int offset = 0; offset < numClasses; offset++) {
		        	num[offset] = Math.exp(num[offset] - max);
		        	denom += num[offset];
		        }
		        Utils.normalize(num, denom);

		        // Update denominator of the gradient of -log(Posterior)
		        double firstTerm;
		        for (int offset = 0; offset < numClasses; offset++) { // Which
		                                                                    // part of x
		        	index = offset * dim;
		        	firstTerm = weights[i] * num[offset];
		        	for (int q = 0; q < dim; q++) {
		        		grad[index + q] += firstTerm * data[i][q];
		        	}
		        }

		        for (int p = 0; p < dim; p++) {
		            grad[cls[i] * dim + p] -= weights[i] * data[i][p];
		        }
		        
			}
		      

			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() > 0){
					for(int pid : n.parents){
						int pOffset = pid*dim;
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						for(int w=0; w<dim; w++){
							grad[nOffset+w] += ADMM_PHO*(sm_x[nOffset+w]-sm_x[pOffset+w]-z[zOffset+w]+u[zOffset+w]);
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=0; w<dim; w++){
						grad[nOffset+w] += ADMM_PHO*(sm_x[nOffset+w]-z[zOffset+w]+u[zOffset+w]);
					}
				}
			}	
		      
			return grad;
		}
		
		
		public double objectiveFunction(double[] currx){
			double nll=0.0;
			int dim = numPredictors+1;
			// update all the internal nodes again
			updateInternalNodes(currx);

			for (int i = 0; i < cls.length; i++) { // ith instance
				double[] exp = new double[numClasses];
				int index;
				for (int offset = 0; offset < numClasses; offset++) {
					index = offset * dim;
					for (int j = 0; j < dim; j++) {
						exp[offset] += data[i][j] * sm_x[index + j];
					}
				}
				double max = exp[Utils.maxIndex(exp)];
		        double denom = 0;
		        double num = exp[cls[i]] - max;
		        
		        for (int offset = 0; offset < numClasses; offset++) {
		        	denom += Math.exp(exp[offset] - max);
		        }

		        nll -= weights[i] * (num - Math.log(denom)); // Weighted NLL
			}
			
			for(Node n : classStructure.allNodes.values()){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() >0){
					for(int pid : n.parents){
						int pOffset = pid*dim;
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						for(int w=0; w<dim; w++){
							nll += (ADMM_PHO/2)*(sm_x[nOffset+w]-sm_x[pOffset+w]-z[zOffset+w]+u[zOffset+w])*(sm_x[nOffset+w]-sm_x[pOffset+w]-z[zOffset+w]+u[zOffset+w]);
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=0; w<dim; w++){
						nll += (ADMM_PHO/2)*(sm_x[nOffset+w]-z[zOffset+w]+u[zOffset+w])*(sm_x[nOffset+w]-z[zOffset+w]+u[zOffset+w]);
					}
				}
			}
			return nll;
		}
		
		

	}

}
