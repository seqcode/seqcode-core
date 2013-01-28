package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.util.ArrayList;
import java.util.Collections;

import Jama.Matrix;

import edu.psu.compbio.seqcode.gse.utils.Pair;


/**
 * DESeqDifferentialEnrichment: implements the differential enrichment method proposed by Anders & Huber (Genome Biology 2010). 
 * i.e. The DESeq method
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class DESeqDifferentialEnrichment extends DifferentialEnrichment{

	protected Normalization normalizer; //DESeq requires a separate normalization method
	protected CountsDataset data;
	protected Lowess[] lowess; //Local regression models. Indexed by condition
	protected double lowessBandwidthPct = (double)0.2;
	//DESeq parameters
	protected double[][] q; 	//Strength parameters (i.e. estimates on condition means). Indexed by unit & condition
	protected double[][] w; 	//Sample variances. Indexed by unit & condition
	protected double[][] z;		//Variance correction parameters. Indexed by unit & condition
	protected double[][] v;		//Variance estimates. Indexed by unit & condition
	protected double[] m; 	//Number of replicates for each condition. Indexed by condition
	
	public DESeqDifferentialEnrichment(Normalization norm){
		super();
		normalizer = norm;
	}
	
	@Override
	public CountsDataset execute(CountsDataset data) {
		this.data = data;
		int ref = data.getFocalCondition();
		//Normalize
		normalizer.normalize(data);
		normalizer.savePairwiseMAPlots(data, true);

		//Initialize DESeq variables
		double[][] k = data.getCounts().getArray();
		double [] scaling = data.scaling;
		int numConds = data.numConds;
		int numUnits = data.numUnits;
		int numSamples = data.numSamples;
		int [] design = data.design;
		m = new double[numConds];
		q = new double[numUnits][numConds];
		w = new double[numUnits][numConds];
		z = new double[numUnits][numConds];
		v = new double[numUnits][numConds];
		lowess = new Lowess[numConds];
		for(int c=0; c<numConds; c++){
			m[c]=0;
			for(int u=0; u<numUnits; u++){
				q[u][c]=0; w[u][c]=0; z[u][c]=0; v[u][c]=0;
			}
		}
		for(int samp=0; samp<numSamples; samp++)
			m[design[samp]]++;
		
		//q: Strength parameters
		for(int u=0; u<numUnits; u++){
			for(int c=0; c<numConds; c++){
				double currQ=0;
				for(int samp=0; samp<numSamples; samp++){
					if(design[samp]==c){
						currQ+=k[u][samp]/scaling[samp];
					}
				}currQ/=m[c];
				q[u][c] = currQ;
			}
		}
		
		//w & z: sample variance estimates
		for(int u=0; u<numUnits; u++){
			for(int c=0; c<numConds; c++){
				if(m[c]>1){
					double currW=0, currInvS=0;
					for(int samp=0; samp<numSamples; samp++){
						if(design[samp]==c){
							currW+=((k[u][samp]/scaling[samp])-q[u][c])*((k[u][samp]/scaling[samp])-q[u][c]);
							currInvS+=1/scaling[samp];
						}
					}
					currW/=(m[c]-1);
					w[u][c] = currW;
					z[u][c] = (q[u][c]/m[c])*currInvS; 
				}
			}
		}

		//Construct graph of (q,w) and estimate smooth local function
		//The function is built here only for conditions that have >1 replicate
		for(int c=0; c<numConds; c++){
			if(m[c]>1){
				ArrayList<QWItem> vals = new ArrayList<QWItem>();
				for(int u=0; u<numUnits; u++)
					vals.add(new QWItem(u, q[u][c], w[u][c], z[u][c]));
				//xval must be monotonically increasing for loess method
				Collections.sort(vals);
				double[] xval = new double[numUnits];
				double[] yval = new double[numUnits];
				Matrix xy = new Matrix(numUnits,2);
				int u=0;
				for(QWItem item : vals){
					xval[u]=(double)item.q;
					yval[u]=(double)item.w;
					xy.set(u,0,item.q);
					xy.set(u,1,item.w);
					u++;
				}
	
				//Run Lowess
				lowess[c]= new Lowess(xval, yval, lowessBandwidthPct);
				double [] fityval = lowess[c].getYEst();
				
				//Variance estimates
				u=0;
				for(QWItem item : vals){
					item.v = fityval[u] - item.z;
					v[item.index][c] = item.v;
					u++;
				}
				
				//Print charts & data
				this.printMeanVarData(xy, data.getCondName(c));
				this.saveMeanVarPlot(xy, fityval, data.getCondName(c), true);
			}
		}

		//If the reference condition has replicates, and other conditions don't, use the reference function for those. 
		if(m[ref]>1){
			for(int c=0; c<numConds; c++){
				if(m[c]==1){
					lowess[c] = lowess[ref];
		}}}

		//If reference condition has only a single replicate, pool all counts and build a single shared reference function. 
		if(m[ref]==1){
			ArrayList<QWItem> vals = new ArrayList<QWItem>();
			for(int c=0; c<numConds; c++)
				for(int u=0; u<numUnits; u++)
					vals.add(new QWItem(u, q[u][c], w[u][c], z[u][c]));
			//xval must be monotonically increasing for loess method
			Collections.sort(vals);
			double[] xval = new double[numUnits];
			double[] yval = new double[numUnits];
			Matrix xy = new Matrix(numUnits,2);
			int u=0;
			for(QWItem item : vals){
				xval[u]=(double)item.q;
				yval[u]=(double)item.w;
				xy.set(u,0,item.q);
				xy.set(u,1,item.w);
				u++;
			}

			//Run Lowess
			for(int c=0; c<numConds; c++){
				lowess[c]= new Lowess(xval, yval, lowessBandwidthPct);
				double [] fityval = lowess[c].getYEst();
			
				//Variance estimates
				u=0;
				for(QWItem item : vals){
					item.v = fityval[u] - item.z;
					v[item.index][c] = item.v;
					u++;
				}
			}
		}
		
		//Update data to include mean & variance estimates
		data.setCondMeans(new Matrix(q));
		data.setCondRawVars(new Matrix(v));
		
		//Calculate p-values for each gene, each condition vs reference
		double[] Kref=new double[numUnits];
		double[] Kx=new double[numUnits]; //Total counts for conditions REF and X
		double[] Qzero=new double[numUnits];   //Pooled mean (under null hypothesis)
		double[] MUref=new double[numUnits];
		double[] MUx=new double[numUnits]; //Mean for conditions REF and X
		double[] SIGMAref; double[] SIGMAx; //Variance for conditions REF and X
		double[] rawSCVref=new double[numUnits]; 
		double[] rawSCVx=new double[numUnits]; //raw SCV
		double[] Pref=new double[numUnits]; double[] Px=new double[numUnits];
		double[] Rref=new double[numUnits]; double[] Rx=new double[numUnits];
		double[] fullVarRef = new double[numUnits]; 
		double[] fullVarX = new double[numUnits];
		
		//Iterate through each condition
		for(int x=0; x<numConds; x++){
			if(x!=ref){
				for(int u=0; u<numUnits; u++){
					//Kref
					Kref[u]=0;
					for(int samp=0; samp<numSamples; samp++){if(design[samp]==ref){Kref[u]+=k[u][samp];}}
			
					//Kx
					Kx[u]=0;
					for(int samp=0; samp<numSamples; samp++){if(design[samp]==x){Kx[u]+=k[u][samp];}}
					
					//Qzero
					Qzero[u]=0;
					for(int samp=0; samp<numSamples; samp++){
						if(design[samp]==x || design[samp]==ref){Qzero[u]+=k[u][samp]/scaling[samp];}
					}Qzero[u]/=(m[ref]+m[x]);
					
					//Means
					MUref[u] = Qzero[u]*(scaling[ref]*m[ref]);
					MUx[u] = Qzero[u]*(scaling[x]*m[x]);
					
					//Set fold
					double fold = q[u][x]>0 ? q[u][ref]/q[u][x] : q[u][ref];
					data.setCondFold(u, x, fold);
					data.setCondMean(u, x, Qzero[u]);
				}
				//Variance
				SIGMAref = lowess[ref].estimateValues(Qzero);
				SIGMAx = lowess[x].estimateValues(Qzero);
				for(int u=0; u<numUnits; u++){
					
					rawSCVref[u] = SIGMAref[u]/(Qzero[u]*Qzero[u]);
					rawSCVx[u] = SIGMAx[u]/(Qzero[u]*Qzero[u]);
						
					fullVarRef[u] = Math.max( MUref[u] + (rawSCVref[u] * Qzero[u] * Qzero[u] * (scaling[ref]*scaling[ref]*m[ref])), MUref[u] * (1+1e-8) );
					fullVarX[u] = Math.max( MUx[u] + (rawSCVx[u] * Qzero[u] * Qzero[u] * (scaling[x]*scaling[x]*m[x])), MUx[u] * (1+1e-8) );
					   
					Pref[u] = MUref[u]/fullVarRef[u];
					Px[u] = MUx[u]/fullVarX[u];
					Rref[u] = (MUref[u]*MUref[u])/(fullVarRef[u]-MUref[u]);
					Rx[u] = (MUx[u]*MUx[u])/(fullVarX[u]-MUx[u]);
					
					if(Pref[u]>0 && Px[u]>0 && Rref[u]>0 && Rx[u]>0){
						Double pObs = 	NegativeBinomialDistrib.pdf( (int)Kref[u], Rref[u], Pref[u] ) * 
										NegativeBinomialDistrib.pdf( (int)Kx[u], Rx[u], Px[u] );
					
						if(pObs.isInfinite() || pObs.isNaN()){
							data.setDEpval(u, x, 1.0);
						}else{
							Pair<Double,Double> totalAbove = addFromMiddle((int)(Kref[u]+Kx[u]), pObs, MUref[u], fullVarRef[u], MUx[u], fullVarX[u], true, 1e-4);
							Pair<Double,Double> totalBelow = addFromMiddle((int)(Kref[u]+Kx[u]), pObs, MUref[u], fullVarRef[u], MUx[u], fullVarX[u], false, 1e-4);
						
							double p = (totalAbove.cdr()+totalBelow.cdr())/(totalAbove.car()+totalBelow.car());
							data.setDEpval(u, x, p);
						}
					}else{
						data.setDEpval(u, x, 1.0);
					}
				}
			}
		}
		
		return data;
	}
	
	/**
	 * Method used to evaluate the upper or lower half of the p(a,b) distribution.
	 * The P value of a pair of observed count sums (kA, kB) is the sum of 
	 * all probabilities less or equal to p(kA, kB), given that the overall sum is kS.
	 * 
	 * This method is ported from pval.c in the DESeq code.
	 *  
	 * @param kS
	 * @param pobs
	 * @param muA
	 * @param vA
	 * @param muB
	 * @param vB
	 * @param upwards
	 * @param eps
	 * @return
	 */
	private Pair<Double,Double> addFromMiddle( int kS, double pobs, double muA, double vA, 
		      double muB, double vB, boolean upwards, double eps){
		int k = (int)((double)kS * muA / (muA+muB));
		if( !upwards )
			k++;
		double total = 0;
		double esttotal = NegativeBinomialDistrib.pdf(kS, (muA+muB)*(muA+muB) / (vA+vB-muA-muB), (muA+muB) / (vA+vB));
		double obstotal = 0;   
		int step = 1;
		double val = 0;
		double sizeA = muA*muA / (vA-muA);
		double probA = muA/vA;
		double sizeB = muB*muB / (vB-muB);
		double probB = muB/vB;
		int knew;
		double prevval;
		while( upwards ? (k < kS) : (k > 0) ) {
			prevval = val;
			while(true){
				if( upwards ) {
					if( k + step > kS )
						step = kS - k;
					knew = k + step; 
				} else {
					if( k - step < 0 )
						step = k;
					knew = k - step; 
				}
				val = NegativeBinomialDistrib.pdf( knew, sizeA, probA) * NegativeBinomialDistrib.pdf( kS - knew, sizeB, probB);
				if( (step == 1) || (Math.abs( val - prevval ) < eps * esttotal/kS)  )
					break;
				step >>= 1;
			}
			k = knew;
			total += val * step;

			if( val <= pobs ) {
				if( prevval <= pobs ) {
					obstotal += val * step;
				} else {
					obstotal += val * ( 1 + (step-1) * ( pobs - val ) / ( prevval - val ) );
				}
			}
			if( (step < Integer.MAX_VALUE) && (Math.abs( val - prevval ) < eps * esttotal/kS / 4) ) {
				step <<= 1;
			}
		}
		return(new Pair<Double,Double>(total, obstotal));
	}

	/**
	 * QWItem: separate class is necessary because the Lowess class requires a sorted list of points, and 
	 * we need a way to map back to the original data structure.  
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	public class QWItem implements Comparable<QWItem>{
		public int index;
		public double q, w;
		public double z, v;
		
		public QWItem(int index, double q, double w, double z){
			this.index=index;
			this.q=q;
			this.w=w;
			this.z=z;
		}
		
		@Override
		public int compareTo(QWItem x) {
			if(this.q<x.q){return -1;}
			else if(this.q>x.q){return 1;}
			return 0;
		}
	}
}
