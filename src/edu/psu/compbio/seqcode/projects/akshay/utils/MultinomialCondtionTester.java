package edu.psu.compbio.seqcode.projects.akshay.utils;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import umontreal.iro.lecuyer.probdistmulti.MultinomialDist;


import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.GFFEntry;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.ParseGFF;

import edu.psu.compbio.seqcode.projects.multigps.framework.BindingManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.MultiGPSConfig;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.sequtils.PointsToEvents;

public class MultinomialCondtionTester {
	
	protected MultiGPSConfig config;
	protected ExperimentManager exptMan;
	protected List<BindingEvent> features;
	protected double minFold = 2;
	protected double genomeLength;
	protected MultinomialDist mnomial;
	
	public MultinomialCondtionTester(MultiGPSConfig con, ExperimentManager em, List<BindingEvent> features, double minFoldChange, double genomeLength ) {
		this.config = con;
		this.exptMan = em;
		this.features = features;
		this.minFold = minFoldChange;
		this.genomeLength = genomeLength;
		
		mnomial = new MultinomialDist(100, new double[2]);
		
		
		
	}
	
	
	public void execute(){
		for(BindingEvent cf : features){
			double[] condSigScaledCounts = new double[exptMan.getConditions().size()];
			
			for(ExperimentCondition ec : exptMan.getConditions()){
				double repSigScaleCount=0;
				for(ControlledExperiment cr : ec.getReplicates()){
					repSigScaleCount = cf.getRepSigHits(cr)/cr.getControlScaling();
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
			String tmpout = "";
			
			for(int i=0; i < condSigScaledCounts.length; i++){
				tmpout = tmpout + condSigScaledCounts[i]+"\t";
			}
			
			System.out.println(tmpout+"\n");
			//System.out.println(cf.getPoint().getLocationString()+"\t"+Double.toString(getMultinomialPval(N,X,P)));
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
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		MultiGPSConfig config = new MultiGPSConfig(gcon, newargs, false);
		econ.setScalingSlidingWindow(50000);
		
		ExperimentManager manager = new ExperimentManager(econ);
		BindingManager bindMan = new BindingManager(manager);
		System.err.println("Conditions:\t"+manager.getConditions().size());
		for(ExperimentCondition c : manager.getConditions()){
			System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
		}
		String siteFile = Args.parseString(args, "pointsFile", null);
		if(siteFile==null){
			 System.exit(1);
		}
		int win = Args.parseInteger(args, "win", 50);
		double minFold = Args.parseDouble(args, "minfold", 2);
		
		List<BindingEvent> features = new ArrayList<BindingEvent>();
		List<Point> Sites = RegionFileUtilities.loadPeaksFromPeakFile(config.getGenome(), siteFile, win);
		 
	
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, bindMan, Sites, win,true);
		features = p2e.execute();
		
		
		MultinomialCondtionTester tester = new MultinomialCondtionTester(config, manager, features, minFold, econ.getMappableGenomeLength());
		tester.execute();
		
		manager.close();
	}
	
	
}
