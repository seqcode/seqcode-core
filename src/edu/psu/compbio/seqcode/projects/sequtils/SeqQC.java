package edu.psu.compbio.seqcode.projects.sequtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariateOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;
import edu.psu.compbio.seqcode.projects.multigps.framework.BackgroundDetector;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;

public class SeqQC {

	Config config;
	Genome gen;
	ExperimentManager manager;
	float[] startcounts =null;
	float densityWindow = 500;
	double testQuantile = 0.1; //Base per-read count distribution on this proportion of the lowest-ranked densities
	int histoMax = 200;
	int poissLowerBound = 1;
	int poissUpperBound = 10; 
	int sesScalingWin=1000;
	boolean verbose=false;
	
	/**
	 * Constructor
	 * @param c
	 * @param man
	 */
	public SeqQC(Config c){
		config = c;
		gen = config.getGenome();
		config.setPerBaseReadFiltering(false);
		config.setSESScaling(true);
		config.setScalingSlidingWindow(sesScalingWin);
		manager = new ExperimentManager(config);
	}
	
	//Accessors & settors
	public void setTestQuantile(double tq){
		if(tq>0 && tq<=1.0)
			testQuantile=tq;
	}
	public void setVerbose(boolean v){verbose=v;}
	
	/**
	 * Run all implemented QC checks
	 */
	public void execute(){
		//Library size calculation
		Map<ControlledExperiment, Double> meanFragmentCoverage = estimateLibrarySize();
		
		//Estimate signal proportion
		Map<ControlledExperiment, Double> sigProps = estimateSignalProportions(meanFragmentCoverage);
	}
	
	/**
	 * estimateLibrarySize
	 * returns the mean coverage of each fragment
	 */
	public Map<ControlledExperiment,Double> estimateLibrarySize(){
		Map<ControlledExperiment, Double> meanFragmentCoverage = new HashMap<ControlledExperiment, Double>();
		
		//If we have multiple experiments, process one at a time
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			System.out.println("Experiment: "+expt.getName()+" = "+String.format("%.1f", expt.getSignal().getHitCount())+ " mapped tags.");
			// 1: Find the density of covered bases surrounding each read
			List<DensityCountPair> densities = new ArrayList<DensityCountPair>();
			Iterator<Region> chroms = new ChromosomeGenerator().execute(gen);
			while (chroms.hasNext()) {
				Region currentRegion = chroms.next();
				//Split the job up into large chunks
	            for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
	                int y = x+config.MAXSECTION; 
	                if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
	                Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
					List<StrandedBaseCount> ipHits = expt.getSignal().getUnstrandedBases(currSubRegion);
					makeHitStartArray(ipHits, currSubRegion, '+');
					float posCounts[] = startcounts.clone();
					makeHitStartArray(ipHits, currSubRegion, '-');
					float negCounts[] = startcounts.clone();
		            
					int halfDWin = (int)densityWindow/2;
					for(int i=currSubRegion.getStart()+halfDWin; i<currSubRegion.getEnd()-halfDWin; i++){  //Note that this means a tiny fraction of reads that are located on the large block boundaries will be ignored
						if(posCounts[i]>0){ //Treat fragments on opposite strands separately
							float dens =0; 
							for(int j=i-halfDWin; j<i+halfDWin; j++){
								if(posCounts[j]>0 && j!=i)
									dens++;
								if(negCounts[j]>0)
									dens++;
							}
							dens /= (densityWindow*2);
							densities.add(new DensityCountPair(dens, posCounts[i]));
						}
						if(negCounts[i]>0){ //Treat fragments on opposite strands separately
							float dens =0; 
							for(int j=i-halfDWin; j<i+halfDWin; j++){
								if(posCounts[j]>0)	
									dens++;
								if(negCounts[j]>0 && j!=i)
									dens++;
							}
							dens /= (densityWindow*2);
							densities.add(new DensityCountPair(dens, negCounts[i]));
						}
					}
	            }
			}
			Collections.sort(densities); //Sort the density pairs in increasing order
			
			//2: Generate a read count per base distribution for the lowest density sites. 
			double currWeight=0, quantileLimit=testQuantile*expt.getSignal().getHitCount();
			RealValuedHistogram histo = new RealValuedHistogram(0,histoMax,histoMax);
			for(DensityCountPair dcp : densities){
				histo.addValue(dcp.getCount());
				currWeight+=dcp.getCount();
				if(currWeight>quantileLimit)
					break;
			}
			
			//3: Fit a distribution to the histogram (Poisson for now)
			DRand re = new DRand();
			int left=poissLowerBound, right=poissUpperBound;
			double xsum=0, xcount=0;
			for(double i=left; i<=right; i++){
				xsum += i*histo.getBin( histo.getBinContainingVal(i));
				xcount += histo.getBin( histo.getBinContainingVal(i));
			}
			double xavg = xsum/xcount;
			UnivariateFunction func = new truncPoisson(xavg, left, right);
			double relativeAccuracy = 1.0e-6;
			double absoluteAccuracy = 1.0e-4;
			UnivariateOptimizer solver = new BrentOptimizer(relativeAccuracy, absoluteAccuracy);
			UnivariatePointValuePair pvp = solver.optimize(100, func, GoalType.MINIMIZE, 0.001, 50.0, xavg);
			double lambda = pvp.getPoint();
			if(verbose)
				System.out.println(String.format("xavg: %.5f\tlambda %.5f", xavg, lambda));
			
			//4: Calculate the library size
			Poisson poiss = new Poisson(lambda, re);
			meanFragmentCoverage.put(expt,  lambda);
			double librarySize = (xcount / (poiss.cdf(right) - poiss.cdf(left - 1))) / testQuantile;
			double libraryCoverage = 1-Math.exp(-expt.getSignal().getHitCount()/librarySize);
			double libraryCoverageDoubledSeq = 1-Math.exp(-(2*expt.getSignal().getHitCount())/librarySize);
			System.out.println("Initial mappable library size = "+String.format("%.1f", librarySize) +" fragments");
			double novelFragsUnderDoubleSeq = (libraryCoverageDoubledSeq-libraryCoverage)*librarySize;
			System.out.println(String.format("Each fragment was sequenced on average %.5f times", lambda));
			System.out.println(String.format("Proportion of library sequenced: %.3f", libraryCoverage));
			System.out.println(String.format("Proportion of library sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", libraryCoverageDoubledSeq, novelFragsUnderDoubleSeq ));
			
			//Print observed and estimated histos
			if(verbose){
				System.out.println("Observed/Expected per-base counts");
				for(int i=0; i<histoMax; i++){
					double obs = histo.getBin(histo.getBinContainingVal(i));
					double exp = librarySize * testQuantile * poiss.pdf(i);
					System.out.println(String.format("%d\t%.0f\t%.0f",i, obs, exp));
				}
			}
		}
		return meanFragmentCoverage;
	}
	
	/**
	 * estimateSignalProportions
	 * @param meanFragmentCoverage
	 */
	public Map<ControlledExperiment,Double> estimateSignalProportions(Map<ControlledExperiment,Double> meanFragmentCoverage){
		Map<ControlledExperiment, Double> sigProps = new HashMap<ControlledExperiment, Double>();
		
		//Estimate proportion via scaling
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			float scaling = meanFragmentCoverage.get(expt).floatValue();
			if(scaling>1)
				expt.correctSignalCounts(scaling);
			sigProps.put(expt, expt.getSigProp());
			System.out.println(String.format("Signal proportion via scaling for: %s = %.5f", expt.getName(), expt.getSigProp()));
		}
		
		//Estimate proportion via distribution-fitting
		BackgroundDetector detector = new BackgroundDetector(config, manager, 200, 100);
		HashMap<Sample, Double> detectedBP = detector.execute();
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			System.out.println(String.format("Signal proportion via fitting for: %s = %.5f", expt.getName(), 1-detectedBP.get(expt.getSignal())));
		}
		
		return sigProps;
	}
	
	//A truncated Poisson function
	private class truncPoisson implements UnivariateFunction {
		protected DRand re = new DRand();
		protected double xavg;
		protected int left, right;
		
		public truncPoisson(double xavg, int left, int right){
			this.xavg = xavg;
			this.left = left;
			this.right = right;
		}
		public double value(double L){
			return -(-Math.log(K(L, left, right, true)) - L + xavg * Math.log(L));		
		}
		public double K(double L, int left, int right, boolean out){
			Poisson poiss = new Poisson(L, re); 
		    if(out && verbose)
		        System.out.println(String.format("K: %.5f %.5f %.5f %d %d",poiss.cdf(right),poiss.cdf(left - 1), L, left, right));
		    return poiss.cdf(right) - poiss.cdf(left - 1);
		}
	}
	
	//Makes integer arrays corresponding to the read starts over the current region
	protected void makeHitStartArray(List<StrandedBaseCount> hits, Region currReg, char strand){
		startcounts = new float[currReg.getWidth()+1];
		for(int i=0; i<=currReg.getWidth(); i++){startcounts[i]=0;}
		for(StrandedBaseCount r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
				startcounts[offset]+=r.getCount();
			}
		}
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	
	/**
	 * Main driver method
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")){
			System.err.println("SeqQC:\n" +
					"\t--geninfo <genome info file>\n" +
					"\t--expt <experiment to test QC>\n" +
					"\t--ctrl <control for this experiment>\n" +
					"\t--format <format of experiment files (SAM/BED/IDX)>\n" +
					"\t--testquantile <proportion of lowest ranked density reads to build distribution from>\n" +
					"\t--verbose\n");
		}else{
			Config con = new Config(args, false);
			
			SeqQC qc = new SeqQC(con);
			if(ap.hasKey("testquantile"))
				qc.setTestQuantile(new Double(ap.getKeyValue("testquantile")));
			if(ap.hasKey("verbose"))
				qc.setVerbose(true);
			qc.execute();
		}
	}
	
	
	private class DensityCountPair implements Comparable<DensityCountPair>{
		private float density;
		private float count;
		
		public DensityCountPair(float d, float c){
			density = d;
			count =c;
		}
		
		public void setCount(float count) {
			this.count = count;
		}
		public void setDensity(float density) {
			this.density = density;
		}
		public float getCount() {
			return count;
		}
		public float getDensity() {
			return density;
		}
		
		// sort according to density
		public int compareTo(DensityCountPair dc) {
			if(density < dc.density)
				return -1;
			else if (density > dc.density)
				return 1;
			else return 0;
		}

	}
}
