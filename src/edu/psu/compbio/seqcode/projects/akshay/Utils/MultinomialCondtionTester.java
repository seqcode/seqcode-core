package edu.psu.compbio.seqcode.projects.akshay.Utils;


import java.util.List;

import umontreal.iro.lecuyer.probdistmulti.MultinomialDist;


import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;

import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;

public class MultinomialCondtionTester {
	
	protected Config config;
	protected ExperimentSet exptSet;
	protected List<BindingEvent> features;
	protected double minFold = 2;
	protected double genomeLength;
	protected MultinomialDist mnomial;
	
	public MultinomialCondtionTester(Config con, ExperimentSet es, List<BindingEvent> features, double minFoldChange, double genomeLength ) {
		this.config = con;
		this.exptSet = es;
		this.features = features;
		this.minFold = minFoldChange;
		this.genomeLength = genomeLength;
		
		mnomial = new MultinomialDist(100, new double[2]);
		
		
		
	}
	
	
	public void execute(){
		for(BindingEvent cf : features){
			double[] condSigScaledCounts = new double[exptSet.getConditions().size()];
			
			for(ExperimentCondition ec : exptSet.getConditions()){
				double repSigScaleCount=0;
				for(ControlledExperiment cr : ec.getReplicates()){
					repSigScaleCount = cr.getSignal().getHitCount()/cr.getControlScaling();
					condSigScaledCounts[ec.getIndex()] += repSigScaleCount;
				}
			}
			
			int N = 0;
			for(int i=0; i<condSigScaledCounts.length; i++){
				N += Math.ceil(condSigScaledCounts[i]);
			}
			
			
			int[] X = new int[condSigScaledCounts.length-1];
			int max_ind = getMaxIndex(condSigScaledCounts);
			
			int count = 0;
			for(int i=0; i< condSigScaledCounts.length; i++){
				if(i!=max_ind){
					X[count] = (int)Math.ceil(condSigScaledCounts[i]);
				}
			}
			
			double[] P = new double[condSigScaledCounts.length-1]; 
			for(int i=0; i<P.length; i++){
				P[i] = 1/(this.minFold + condSigScaledCounts.length-1);
			}
			
			System.out.println(cf.getPoint().getLocationString()+"\t"+Double.toString(getMultinomialPval(N,X,P)));
		}
	}
	
	public double getMultinomialPval(int n, int[]X, double[] P){
		double ret =0.0;
		this.mnomial.setParams(n, P);
		ret = mnomial.cdf(X);
		return ret;
	}
	
	public int getMaxIndex(double[] list){
		int ret=0;
		double maax = Double.MIN_VALUE;
		for(int i=0; i<list.length; i++){
			if(list[i] > maax){
				maax = list[i];
				ret = i;
			}
		}
		return ret;
	}
	
	
	public static void main(String[] args){
		//Hack to set a default fixedpb limit. Need a better way to set per-application defaults
		String [] newargs = new String[args.length+2];
		for(int a=0; a<args.length; a++)
			newargs[a] = args[a];
		newargs[newargs.length-2] = "--fixedpb";
		newargs[newargs.length-1] = "1000";
		
		Config config = new Config(newargs, false);
		config.setMedianScaling(true);
		config.setScalingSlidingWindow(50000);
		
		ExperimentManager manager = new ExperimentManager(config);
		
		
	}
	
	
}
