package edu.psu.compbio.seqcode.projects.akshay.MultiSeq;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import weka.core.Optimization;
import weka.core.RevisionUtils;
import weka.core.Utils;
import edu.psu.compbio.seqcode.projects.akshay.MultiSeq.SeqUnwinder.ClassRelationStructure;
import edu.psu.compbio.seqcode.projects.akshay.MultiSeq.SeqUnwinder.ClassRelationStructure.Node;

/** 
 * The main class that does all the training for SeqUnwinder. SeqUnwinder implements an L1 norm Logistic regression model with multiple classes related by a given class 
 * relationships. 
 * This class learns the SeqUnwinder model using the ADMM framework. 
 * The x-update subproblem in the ADMM framework is done using the BGFS method.
 * 
 * @author akshaykakumanu
 *
 */

public class LOne extends Optimizer {
	
	// Fixed ADMM parameters
	
	/** Relaxation parameter (to help faster convergence) */
	public final double ADMM_ALPHA = 1.9;
	/** Absolute feasibility tolerance for the primal and dual feasibility conditions */
	public final double ADMM_ABSTOL = 1E-4; 
	/** Relative  feasibility tolerance for the primal and dual feasibility conditions */
	public final double ADMM_RELTOL = 1E-2;
	
	// Fixed SeqUnwinder parameters
	
	/** Tolerence for internal Nodes convergence */
	public final double NODES_tol = 1E-4;
	
	// Tunable ADMM parameters
	
	/** The maximum number of allowed iterations for the ADMM algorithm */
	public int ADMM_maxItr = 30;
	/** Augmented Lagrangian parameter rho */
	public double ADMM_pho = 0.001; 
	/** The maximum allowed value for pho */
	public double ADMM_pho_max = 100000000;
	
	// BGFS parameters 
	
	/** The maximum number of allowed iterations for the BGFS algorithm */
	public int BGFS_maxIts=-1;
	
	// SeqUnwinder parameters
	/** The total number of predictors/features in the model (does not include the intercept term) */
	protected int numPredictors;
	/** Total number of classes to be predicted */
	protected int numClasses;
	/** Relationships between the different nodes in SeqUnwinder */
	protected ClassRelationStructure classStructure;
	/** Total number of nodes (internal nodes and classes) in SeqUnwinder */
	protected int numNodes;
	/** The L1 regularization parameter in the SeqUnwinder model */
	protected double regularization;
	/** Maximum number of iterations to update internal nodes. For small number of node levels (usually 3 to 4), we might need very few iterations*/
	protected int NODES_maxItr=10;
	/** Internal flag that indicated if ADMM has been run before */
	protected boolean ranADMM=false;
	
	// SeqUnwinder and ADMM variables
	
	/** Current feature weights for all the nodes in the SeqUnwinder (also contains the intercept term). Dimension :- (numPredictors+1)*numNodes */
	public double[] sm_x;
	/** Current feature weights for all the leaf nodes (classes) in the SeqUnwinder (also contains the intercept term). Dimension :- (numPredictors+1)*numClasses */
	public double[] x; 
	/** Current values of z (z-step in ADMM). Dimension :- (numPredictors+1)*numNodes*numNodes */
	public double[] z; // for the Z-step
	/** Value of z at previous iteration. Needed to assess convergence  Dimension :- (numPredictors+1)*numNodes*numNodes */
	public double[] zold;
	/** Current values of the augmented lagrange dual variables. Dimension :- (numPredictors+1)*numNodes*numNodes  */
	public double[] u;
	/** Stores the primal residuals over the course of the ADMM algorithm Dimension:- [ADMM_maxItr][numNodes*numNodes] */
	public double[][] history_primal;
	/** Stores the dual residuals over the course of the ADMM algorithm Dimension:- [ADMM_maxItr][numNodes*numNodes] */
	public double[][] history_dual;
	/** Boundry conditions needed for the the newton's method. However, null in this case */
	public double[][] b;
	
	
	// SeqUnwinder training data
	
	/** Training data (Instances)*/
	public double[][] data;
	/** Weights of the instances */
	protected double[] weights;
	/** Instance class membership */
	protected int[] cls;
	
	// Misc
	
	/** Boolean variable indicating debug mode */
	protected boolean sm_Debug;
	
	
	//Settors
	public void setBGFSmaxItrs(int m){BGFS_maxIts=m;}
	public void setADMMmaxItrs(int m){ADMM_maxItr = m;}
	public void setSeqUnwinderMaxIts(int m){NODES_maxItr = m;}
	public void setInstanceWeights(double[] w){weights=w;}
	public void setClsMembership(int[] c){cls=c;}
	public void setNumPredictors(int p){numPredictors=p;}
	public void setNumClasses(int c){numClasses = c;}
	public void setClassStructure(ClassRelationStructure rel){classStructure = rel; setNumNodes(rel.allNodes.size());}
	public void setNumNodes(int n){numNodes = n;}
	public void setRidge(double r){regularization = r;}
	public void setDebugMode(boolean debug){sm_Debug =debug;}
	public void setPho(double ph){ADMM_pho = ph;}
	
	//gettors
	public double[] getX(){return x;}
	public double[] getZ(){return z;}
	public double[] getU(){return u;}
	public double[] getsmX(){return sm_x;}
	
	
	// Initialize
	public void initZandU(){
		int dim = numPredictors+1;
		z= new double[numNodes*numNodes*dim];
		zold = new double[numNodes*numNodes*dim];
		u= new double[numNodes*numNodes*dim];
		history_primal = new double[ADMM_maxItr][numNodes*numNodes];
		history_dual = new double[ADMM_maxItr][numNodes*numNodes];
	}
	
	public LOne(double[] xinit, double[] sm_xinit, double[][] d, double[][] bc) {
		x = xinit;
		sm_x = sm_xinit;
		data=d;
		b=bc;
	}
	
	public void execute() throws Exception{
		
		for(int it=0; it<NODES_maxItr; it++){
			
			double[] sm_x_old = new double[sm_x.length];
			
			for(int i=0; i<sm_x.length; i++){
				sm_x_old[i] = sm_x[i];
			}
			
			// First, run admm on leaf nodes
			executeADMM();
			
			// Now update the internal nodes
			updateInternalNodes();
			
			// Check Convergence
			boolean converged = true;
			for(Node n : classStructure.allNodes.values()){
				double diff = 0.0;
				for(int w=0; w<(numPredictors+1); w++){
					diff += Math.pow(sm_x_old[n.nodeIndex*(numPredictors+1)+w]-sm_x[n.nodeIndex*(numPredictors+1)+w],2);
				}
				diff = Math.sqrt(diff);
				if( diff > Math.sqrt(numPredictors)*NODES_tol){
					converged=false;
					break;
				}
			}
			
			if(converged)
				break;
		}
		
	}
	
	public void executeADMM() throws Exception{
		int dim = numPredictors+1;
		for(int itr=0; itr<ADMM_maxItr; itr++){
			System.err.print(". "+ itr + " .");
			
			
			// Update pho 
			if(itr >0 && !ranADMM && ADMM_pho < ADMM_pho_max)
				updatePhoAndU(itr-1);
			
			System.err.print(" "+ADMM_pho+" ");
			
			// Update x
			//double xmin = updateX();
			updateApproxX();
			
			// Copy z to zold
			for(int i=0; i<z.length; i++){
				zold[i] = z[i];
			}
			
			// Calculate over-relaxed x:- xrel
			double[] xhat = new double[z.length];
			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size()>0){
					for(int pid : n.parents){
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						int pOffset = pid*dim;
						for(int w=0; w<dim; w++){
							xhat[zOffset+w] = ADMM_ALPHA*(sm_x[nOffset+w])+(1-ADMM_ALPHA)*zold[zOffset+w]-ADMM_ALPHA*sm_x[pOffset+w];
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=0; w<dim; w++){
						//xrel[zOffset+w] = ADMM_ALPHA*sm_x[nOffset+w]+(1-ADMM_ALPHA)*zold[zOffset+w]+u[zOffset+w];
						xhat[zOffset+w] = ADMM_ALPHA*sm_x[nOffset+w]+(1-ADMM_ALPHA)*zold[zOffset+w];
					}
				}
			}
			
			for(int i=0; i<z.length;i++){
				z[i] = xhat[i]+u[i];
			}
			
			// Z-update
			double zmin = updateZ((2*regularization)/ADMM_pho);
			
			//U-update
			
			for(int i=0; i<u.length; i++){
				u[i] = u[i] + xhat[i] - z[i]; 
			}
			
			
			// Check Convergence
			boolean converged = hasADMMConverged(itr);
			
			// Print the prinal and dual residuals
			
			if(sm_Debug){
				double primal = 0.0;
				double dual = 0.0;
				for(int i=0; i<(numNodes*numNodes); i++){
					primal += history_primal[itr][i];
					dual += history_dual[itr][i];
				}
				System.err.println("Primal residual "+ primal + " , Dual residual "+ dual);
			}
						
			
			if(converged){
				System.err.println();
				System.err.println("ADMM has converged after "+itr+" iterations !!");
				ranADMM=true;
				break;
			}
		}
		ranADMM = true;
		
	}
	
	
	public void updatePhoAndU(int A_itr){
		
		double primal_residuals=0; // Sum of all the primal residuals
		double dual_residuals=0; // Sum of all the dual residuals
		
		double mu = 10; // maintains the primal and dual residuals within a factor of mu of one another
		double tao=2; // the factor by with pho will be increased or decreased at each iterations
		
		for(int i=0; i<(numNodes*numNodes); i++){
			primal_residuals += history_primal[A_itr][i];
			dual_residuals += history_dual[A_itr][i];
		}
		
		if(primal_residuals > mu*dual_residuals){ // if primal residual are greater than dual by a factor of mu; decrease pho 
			double old_pho = ADMM_pho;
			ADMM_pho = Math.min(ADMM_pho_max,ADMM_pho*tao);
			double fold =  ADMM_pho/old_pho;
			for(int i=0; i<u.length; i++){ // update u
				u[i] = u[i]/fold;
			}
		}else if(dual_residuals > primal_residuals*mu){
			ADMM_pho = ADMM_pho/tao;
			for(int i=0; i<u.length; i++){ // update u
				u[i] = u[i]*tao;
			}
		}
		
	}
	
	public boolean hasADMMConverged(int itr){
		boolean converged = true;
		int dim = numPredictors + 1;
		
		double[] primal_tol = new double[numNodes*numNodes];
		double[] dual_tol = new double[numNodes*numNodes];
		
		for(Node n : classStructure.leafs){
			double xnorm = getL2NormX(n.nodeIndex);
			if(n.parents.size() > 0){ // If this node has parents
				for(int pid : n.parents){ // Over all the parents of this node
					double znorm = getL2NormZ(n.nodeIndex, pid);
					double unorm = getL2NormU(n.nodeIndex, pid);
					double cnorm = getL2NormX(pid);
					primal_tol[n.nodeIndex*numNodes+pid] = Math.sqrt(dim)*ADMM_ABSTOL + ADMM_RELTOL*Math.max(xnorm, Math.max(znorm, cnorm));
					dual_tol[n.nodeIndex*numNodes+pid] = Math.sqrt(dim)*ADMM_ABSTOL + ADMM_RELTOL*Math.sqrt(ADMM_pho)*unorm;
					
					double deltaZnorm = 0;
					double primalResidualNorm = 0;
					int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
					for(int w=0; w<dim; w++){
						deltaZnorm += Math.pow(ADMM_pho*(z[zOffset+w] - zold[zOffset+w]), 2);
						primalResidualNorm += Math.pow((sm_x[n.nodeIndex*dim+w]-z[zOffset+w]-sm_x[pid*dim+w]),2);
					}
					deltaZnorm = Math.sqrt(deltaZnorm);
					primalResidualNorm = Math.sqrt(primalResidualNorm);
					history_primal[itr][n.nodeIndex*numNodes+pid] = primalResidualNorm;
					history_dual[itr][n.nodeIndex*numNodes+pid] = deltaZnorm;
				}
			}else{
				double znorm = getL2NormZ(n.nodeIndex, n.nodeIndex);
				double unorm = getL2NormU(n.nodeIndex, n.nodeIndex);
				primal_tol[n.nodeIndex*numNodes+n.nodeIndex] = Math.sqrt(dim)*ADMM_ABSTOL + ADMM_RELTOL*Math.max(xnorm, znorm);
				dual_tol[n.nodeIndex*numNodes+n.nodeIndex] = Math.sqrt(dim)*ADMM_ABSTOL + ADMM_RELTOL*Math.sqrt(ADMM_pho)*unorm;
				double deltaZnorm = 0;
				double primalResidualNorm = 0;
				int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
				for(int w=0; w<dim; w++){
					deltaZnorm += Math.pow(ADMM_pho*(z[zOffset+w] - zold[zOffset+w]), 2);
					primalResidualNorm += Math.pow((sm_x[n.nodeIndex*dim+w]-z[zOffset+w]),2);
				}
				deltaZnorm = Math.sqrt(deltaZnorm);
				primalResidualNorm = Math.sqrt(primalResidualNorm);
				history_primal[itr][n.nodeIndex*numNodes+n.nodeIndex] = primalResidualNorm;
				history_dual[itr][n.nodeIndex*numNodes+n.nodeIndex] = deltaZnorm;
			}
		}
		
		for(int i=0; i<primal_tol.length; i++){
			if(history_primal[itr][i] > primal_tol[i])
				converged=false;
			if(history_dual[itr][i] > dual_tol[i])
				converged=false;
			if(!converged)
				break;
		}
		
		return converged;
	}
	
	
	
	public double updateZ(double pho){
		
		int dim = numPredictors+1;
		
		double z_ret=0;
		
		for(int i=0; i<z.length; i++){
			z[i] = z[i] - Math.signum(z[i])*Math.min(pho, Math.abs(z[i]));
			//z[i] = Math.max(0, xrel[i]-pho) - Math.max(0, -xrel[i] - pho);
		}
		
		// Compute the found minimum value for the z-update step
		// Sum over all the Znp update values
		
		/*
		for(Node n : classStructure.leafs){
			double znp=0;
			if(n.parents.size() > 0){ // if the leaf node has parents
				for(int pid : n.parents){ // over all the parents of the current node
					// First part, one-norm part
					double firstPart = 0;
					int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
					for(int w=0; w<dim; w++){
						firstPart += regularization*Math.abs(z[zOffset+w]);
					}
					double secondPart = 0;
					int nOffset = n.nodeIndex*dim;
					int pOffset = pid*dim;
					for(int w=0; w<dim; w++){
						secondPart += (ADMM_pho/2)*Math.pow(sm_x[nOffset+w]-z[zOffset+w]-sm_x[pOffset+w]+u[zOffset+w], 2);
					}
					znp = firstPart + secondPart;
				}
			}else{ // if the lead node has no parents
				double firstPart = 0;
				int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
				int nOffset = n.nodeIndex*dim;
				for(int w=0; w<dim; w++){
					firstPart += regularization+Math.abs(z[zOffset+w]);
				}
				double secondPart = 0;
				for(int w=0; w<dim; w++){
					secondPart += (ADMM_pho/2)*Math.pow(sm_x[nOffset+w]-z[zOffset+w]+u[zOffset+w], 2);
				}
				znp = firstPart + secondPart;
			}
			z_ret += znp;
		}
		*/
		return z_ret;
	}
	
	/**
	 * 
	 * @return the current value of the function that is minimized during the x-update 
	 * @throws Exception
	 */
	public double updateX() throws Exception{
		
		double nll_ret = 0;
		
		// First, create and opt object for the BGFS algorithm
		OptObject oO = new OptObject();
		Optimization opt = new OptEng(oO);
		
		if(BGFS_maxIts == -1){ // Search until convergence
			x = opt.findArgmin(x, b);
			while (x == null) {
				x = opt.getVarbValues();
				x = opt.findArgmin(x, b);
			}
			nll_ret = opt.getMinFunction();
		}else{
			opt.setMaxIteration(BGFS_maxIts);
			x = opt.findArgmin(x, b);
			if (x == null) {
				x = opt.getVarbValues();
			}
			nll_ret = opt.getMinFunction();
		}
			
		// Now copy the leaf node weights (i.e x) to sm_x
		for(Node n: classStructure.leafs){
			int dim = numPredictors+1;
			int nOffset = n.nodeIndex*dim;
			for(int w=0; w<dim; w++){
				sm_x[nOffset+w] = x[nOffset+w];
			}
		} 
		return nll_ret;
	}
	
	public void updateApproxX() throws Exception {
		
		OptObject oO = new OptObject();
		
		int[] iflag = new int[1];
		double obj= oO.objectiveFunction(x);
		double[] grad = oO.evaluateGradient(x);
		int m = 5;
		double[] diag = new double[x.length];
		int[] iprint = new int[2];
		double eps = 0.001;
		double xtol = 10e-16;
		
		LBFGS.lbfgs(x.length, m, x, obj, grad, false, diag, iprint, eps, xtol, iflag);
		
		while(iflag[0] == 1 ){
			//re-evaluate the objective and the gradient
			obj = oO.objectiveFunction(x);
			grad = oO.evaluateGradient(x);
			
			LBFGS.lbfgs(x.length, m, x, obj, grad, false, diag, iprint, eps, xtol, iflag);
		}
		
		// Now copy the leaf node weights (i.e x) to sm_x
		for(Node n: classStructure.leafs){
			int dim = numPredictors+1;
			int nOffset = n.nodeIndex*dim;
			for(int w=0; w<dim; w++){
				sm_x[nOffset+w] = x[nOffset+w];
			}
		}
	
	}
	
	
	
	// Slave methods
	
	private void updateInternalNodes(){
		// First update odd layaer nodes
		for(int l=1; l<classStructure.numLayers; l+=2){
			for(Node n : classStructure.layers.get(l)){// Get nodes in this layer
					updateNode(n);
			}
		}
		// Now update even layer nodes except the leaf node
		for(int l=2; l<classStructure.numLayers; l+=2){
			for(Node n : classStructure.layers.get(l)){// Get nodes in this layer
					updateNode(n);
			}
		}
	}
	
	private void updateNode(Node n){
		int dim = numPredictors+1;
		int nOffset = n.nodeIndex*dim;
		// Note the intercept term not included
		for(int w=1; w<dim; w++){
			List<Double> xs = new ArrayList<Double>();
			for(int pid : n.parents){
				xs.add(sm_x[pid*dim+w]);
			}
			for(int cid : n.children){
				xs.add(sm_x[cid*dim+w]);
			}
			Collections.sort(xs);
			int midInd = xs.size()/2;
			sm_x[nOffset+w] = xs.get(midInd);
		}
	}
	
	private double getL2NormX(int nodeIndex){
		double norm=0;
		int dim = numPredictors+1;
		int offset = nodeIndex*dim;
		for(int w=0; w<dim; w++){
			norm += Math.pow(sm_x[offset+w], 2);
		}
		return Math.sqrt(norm);
	}
	
	private double getL2NormZ(int nInd, int pInd){
		double norm=0;
		int dim = numPredictors+1;
		int zOffset = (nInd*numNodes*dim)+(pInd*dim);
		for(int w=0; w<dim; w++){
			norm += Math.pow(z[zOffset+w], 2);
		}
		return Math.sqrt(norm);
	}
	
	private double getL2NormU(int nInd, int pInd){
		double norm=0;
		int dim = numPredictors+1;
		int uOffset = (nInd*numNodes*dim)+(pInd*dim);
		for(int w=0; w<dim; w++){
			norm += Math.pow(u[uOffset+w], 2);
		}
		return Math.sqrt(norm);
	}
	
	
	/**
	 * This class implements two things:-
	 * It calculates the gradient for the x-update sub-problem. (The BGFS method will need this)
	 * It calculates the overall objective function for the x-update subproblem. (The BGFS method will need this)
	 * @author akshaykakumanu
	 *
	 */
	public class OptObject {
		
		/**
		 * Claclulates the gradient
		 * @param currx
		 * @return
		 */
		public double[] evaluateGradient(double[] x){
			
			double[] grad = new double[x.length];
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
		        		exp += data[i][j]*x[index + j];
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
						for(int w=1; w<dim; w++){
							grad[nOffset+w] += ADMM_pho*(x[nOffset+w]-sm_x[pOffset+w]-z[zOffset+w]+u[zOffset+w]);
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=1; w<dim; w++){
						grad[nOffset+w] += ADMM_pho*(x[nOffset+w]-z[zOffset+w]+u[zOffset+w]);
					}
				}
			}	
		      
			if(sm_Debug){
				//System.err.println(grad[20]);
				//System.err.println(grad[dim+20]);
				//System.err.println(grad[2*dim+20]);
			}
			
			return grad;
		}
		
		/**
		 * Claclulates the objective function
		 * @param currx
		 * @return
		 */
		public double objectiveFunction(double[] x){
			double nll=0.0;
			int dim = numPredictors+1;

			for (int i = 0; i < cls.length; i++) { // ith instance
				double[] exp = new double[numClasses];
				int index;
				for (int offset = 0; offset < numClasses; offset++) {
					index = offset * dim;
					for (int j = 0; j < dim; j++) {
						exp[offset] += data[i][j] * x[index + j];
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
			
			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() >0){
					for(int pid : n.parents){
						int pOffset = pid*dim;
						int zOffset = (n.nodeIndex*numNodes*dim)+(pid*dim);
						for(int w=1; w<dim; w++){
							nll += (ADMM_pho/2)*(x[nOffset+w]-sm_x[pOffset+w]-z[zOffset+w]+u[zOffset+w])*(x[nOffset+w]-sm_x[pOffset+w]-z[zOffset+w]+u[zOffset+w]);
						}
					}
				}else{
					int zOffset = (n.nodeIndex*numNodes*dim)+(n.nodeIndex*dim);
					for(int w=1; w<dim; w++){
						nll += (ADMM_pho/2)*(x[nOffset+w]-z[zOffset+w]+u[zOffset+w])*(x[nOffset+w]-z[zOffset+w]+u[zOffset+w]);
					}
				}
			}
			
			if(sm_Debug){
				//System.err.println("Negative Log Likelihood: "+nll);
			}
			
			return nll;
		}
		
	}
	
	public class OptEng extends Optimization {
		
		OptObject m_oO = null;
		
		public OptEng(OptObject oO){
			m_oO = oO;
		}
		 
		@Override
		public String getRevision() {
			return RevisionUtils.extract("$Revision: 11247 $");
		}

		@Override
		protected double objectiveFunction(double[] x) throws Exception {
			return m_oO.objectiveFunction(x);
		}

		@Override
		protected double[] evaluateGradient(double[] x) throws Exception {
			return m_oO.evaluateGradient(x);
		}
	}
	
	// To clear memory
	public void clearOptimizer(){
		data=null;
		sm_x=null;
		x=null;
		z=null;
		u=null;
	}
}
