package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

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
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;

/**
 * BackgroundDetector:   
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class BackgroundDetector {

	protected ExperimentManager manager; 
	protected Config config;
	protected Genome gen;
	protected float binWidth=0, binStep, winExt;
	protected boolean stranded=false;
	protected HashMap<Sample, RealValuedHistogram> sampleHistos =new HashMap<Sample, RealValuedHistogram>();
	protected double[] sampleTotals;
	protected int histoMax = 500;
	protected int poissUpperBound = 50;
	protected int poissUpperBoundMin = 10;
	protected double cdfPercOfUniform = 0.9; //Percentage of uniform-assumption CDF - used to set upper bound on truncated Poisson
	protected List<Region> regionsToIgnore;
	
	public BackgroundDetector(Config c, ExperimentManager man, float binW, float binS){
		manager = man;
		config = c; 
		gen = config.getGenome();
		binWidth = binW;
		binStep = binS;
		winExt = binWidth/2;
		
		regionsToIgnore = config.getRegionsToIgnore();
		
		sampleTotals =new double[manager.getExperimentSet().getSamples().size()];
		for(int s=0; s<sampleTotals.length; s++)
    		sampleTotals[s]=0;
		
		//Initialize histograms
		for(Sample samp : manager.getExperimentSet().getSamples())
    		sampleHistos.put(samp, new RealValuedHistogram(0,histoMax,histoMax));
    	
		for(Sample s : manager.getExperimentSet().getSamples())
			if(s!=null){
				System.err.println("Sample "+s.getName()+"\t"+s.getHitCount());
				System.err.println("Mean if uniform: "+s.getHitCount()/(gen.getGenomeLength()/binWidth));
			}
	}
	
	

	/**
	 * Calculate binned coverage histograms for each sample and fit truncated Poissons to the lower end of the histograms
	 * Returns a hash map of samples to background proportions. 
	 */
	public HashMap<Sample, Double> execute(){
		Iterator<Region> testRegions = new ChromosomeGenerator().execute(config.getGenome());
		
		Thread[] threads = new Thread[config.getMaxThreads()];
        ArrayList<Region> threadRegions[] = new ArrayList[config.getMaxThreads()];
        int i = 0;
        for (i = 0 ; i < threads.length; i++) {
            threadRegions[i] = new ArrayList<Region>();
        }
        while(testRegions.hasNext()){
        	Region r = testRegions.next(); 
            threadRegions[(i++) % config.getMaxThreads()].add(r);
        }

        for (i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new PotentialRegionFinderThread(threadRegions[i]));
            t.start();
            threads[i] = t;
        }
        boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) { }
            for (i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }
        }
        
        //Fit the Poissons
        HashMap<Sample, Double> backProps = new HashMap<Sample, Double>();
        for(Sample samp : manager.getExperimentSet().getSamples()){
        	if(samp!=null){
        		backProps.put(samp, fitPoisson(sampleHistos.get(samp), samp));
        	}
        }
        return backProps;
	}
	
	/**
	 * Print the contents of the histograms
	 */
	public void print(){
		for(Sample samp : manager.getExperimentSet().getSamples()){
			if(samp!=null){
				System.out.println(samp.getName());
				System.out.println(sampleHistos.get(samp).contentsToString());
			}
		}
	}
	
	/**
	 * Fit a truncated Poisson to the contents of a histogram.
	 * Returns the background proportion in the sample
	 *
	 */
	public double fitPoisson(RealValuedHistogram h, Sample samp){
		DRand re = new DRand();
		//Heuristic to find the upper bound for the truncated Poisson
		int pUpper = poissUpperBound;
		double uniformMean = samp.getHitCount()/(gen.getGenomeLength()/binWidth);
		if(cdfPercOfUniform>0 && cdfPercOfUniform<=1){
			Poisson uniPoiss = new Poisson(uniformMean, re);
			double tmpProp=0;
			int i=0;
			while(tmpProp<cdfPercOfUniform){
				tmpProp = uniPoiss.cdf(i);
				i++;
			}
			pUpper=Math.max(i,poissUpperBoundMin);
		}
		System.out.println("Truncated Poisson Upper Bound:\t"+pUpper);
		
		//Fit the Poisson
		int left=0, right=pUpper;
		double xsum=0, xcount=0;
		for(double i=left; i<=right; i++){
			xsum += i*h.getBin( h.getBinContainingVal(i));
			xcount += h.getBin( h.getBinContainingVal(i));
		}
		double xavg = xsum/xcount;
		UnivariateFunction func = new truncPoisson(xavg, left, right);
		double relativeAccuracy = 1.0e-6;
		double absoluteAccuracy = 1.0e-4;
		UnivariateOptimizer solver = new BrentOptimizer(relativeAccuracy, absoluteAccuracy);
		UnivariatePointValuePair pvp = solver.optimize(100, func, GoalType.MINIMIZE, 0.001, 50.0, xavg);
		double lambda = pvp.getPoint();
		System.out.println("xavg: "+ xavg+"\tlambda: "+lambda);
		
		//Calculate the background proportion
		Poisson poiss = new Poisson(lambda, re);
		double backsize = xsum / (poiss.cdf(right) - poiss.cdf(left - 1));
		double backprop = Math.min(backsize / sampleTotals[samp.getIndex()], 1.0);
		System.out.println("Background= "+ backsize+" / "+sampleTotals[samp.getIndex()]+" =\t"+backprop);
		
		return backprop;
	}
	
	
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
		    if(out)
		        System.out.println("K: "+poiss.cdf(right)+" "+poiss.cdf(left - 1)+" "+ L+" "+ left+" "+ right);
		    return poiss.cdf(right) - poiss.cdf(left - 1);
		}
	}
		
    class PotentialRegionFinderThread implements Runnable {
        private Collection<Region> regions;
        private double[][] starts=null;
        private List<Region> threadPotentials = new ArrayList<Region>();
        
        public PotentialRegionFinderThread(Collection<Region> r) {
            regions = r;
        }
        
        public void run() {
        	HashMap<Sample, RealValuedHistogram> tmpSampleHistos =new HashMap<Sample, RealValuedHistogram>();
        	double[] tmpSampleTotals =new double[manager.getExperimentSet().getSamples().size()];
        	for(Sample samp : manager.getExperimentSet().getSamples())
        		if(samp!=null){
        			tmpSampleHistos.put(samp, new RealValuedHistogram(0,histoMax,histoMax));
        			tmpSampleTotals[samp.getIndex()]=0;
        		}
        	
        	int expansion = (int)winExt;
        	for (Region currentRegion : regions) {
            	Region lastPotential=null;
                //Split the job up into large chunks
                for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
                    int y = (int) (x+config.MAXSECTION+(expansion)); //Leave a little overhang to handle enriched regions that may hit the border. Since lastPotential is defined above, a region on the boundary should get merged in.
                    if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                    Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
                    
                    List<List<StrandedBaseCount>> hits = new ArrayList<List<StrandedBaseCount>>();
                    
                    synchronized(manager){
	                    //Initialize the read lists
                    	for(Sample samp : manager.getExperimentSet().getSamples()){
                    		hits.add(new ArrayList<StrandedBaseCount>());
                    	}
                    	//Load reads by replicate
                    	for(Sample samp : manager.getExperimentSet().getSamples()){
                    		if(samp!=null)
                    			hits.get(samp.getIndex()).addAll(samp.getUnstrandedBases(currSubRegion));
                    	}
                    }
            		int numStrandIter = stranded ? 2 : 1;
                    for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
                        //If stranded peak-finding, run over both strands separately
                        char str = !stranded ? '.' : (stranditer==1 ? '+' : '-');
					 
                        makeStartLandscape(hits, currSubRegion, binWidth, binStep, str);
                        double binnedStarts[][] = starts.clone();
					
                        //Scan regions
                        int currBin=0;
                        for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)binWidth; i+=(int)binStep){
                        	boolean regionPasses=false;
                        	for(Sample samp : manager.getExperimentSet().getSamples()){
                        		if(samp!=null){
                        			double winHits=binnedStarts[samp.getIndex()][currBin];
                        			tmpSampleHistos.get(samp).addValue(winHits);
                        			tmpSampleTotals[samp.getIndex()]+=winHits;
                        		}
                        	}
                            currBin++;
                        }
					}
                }
            }
        	synchronized(sampleHistos){
        		for(Sample samp : manager.getExperimentSet().getSamples()){
        			if(samp!=null)
        				sampleHistos.get(samp).addHistogram(tmpSampleHistos.get(samp));
        		}
        	}	
        	synchronized(sampleTotals){
        		for(Sample samp : manager.getExperimentSet().getSamples()){
        			if(samp!=null)
        				sampleTotals[samp.getIndex()]+=tmpSampleTotals[samp.getIndex()];
        		}
        	}
        }
        

		//Makes integer array corresponding to the binned read start landscape over the current region.
        //No needlefiltering here as that is taken care of during read loading (i.e. in Sample)
    	protected void makeStartLandscape(List<List<StrandedBaseCount>> hits, Region currReg, float binWidth, float binStep, char strand){
    		int numBins = (int)(currReg.getWidth()/binStep);
    		starts = new double[hits.size()][numBins+1];
    		for(Sample samp : manager.getExperimentSet().getSamples()){
    			if(samp!=null){
	            	List<StrandedBaseCount> currHits = hits.get(samp.getIndex());
		    		for(int i=0; i<=numBins; i++){ starts[samp.getIndex()][i]=0; }
		    		for(StrandedBaseCount r : currHits){
		    			if(strand=='.' || r.getStrand()==strand){
		    				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
		    				int binstart = inBounds((int)((double)offset/binStep)-1, 0, numBins);
		    				int binend = inBounds((int)((double)offset/binStep), 0, numBins);
		    				for(int b=binstart; b<=binend; b++)
		    					starts[samp.getIndex()][b]+=r.getCount();
		    			}
		    		}
    			}
        	}
    	}
    	protected final int inBounds(int x, int min, int max){
    		if(x<min){return min;}
    		if(x>max){return max;}
    		return x;
    	}
    	
    }
    
    
    /**
	 * This main method is for testing the BackgroundDetector
	 * @param args
	 */
	public static void main(String[] args){
		
		Config config = new Config(args, false);
		
		if(config.helpWanted()){
			System.err.println("BackgroundDetector:");
			System.err.println("Genome:" +
					"\t--species <Organism;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiments:\n" +
					"\t--expt <read file name> AND --format <SAM/BED/IDX/BOWTIE/NOVO/ELAND>\n" +
					"AND/OR" +
					"\t--rdbexpt <ReadDB experiment identifier>\n" +
					"Other:\n" +
					"\t--threads <number of threads to use>\n" +
					"\t--binwidth <bin width>\n" +
					"\t--binstep <bin step>\n" +
					"\t--fixedpb <fixed per base limit>\n" +
					"\t--poissongausspb <filter per base using Poisson Gaussian sliding window>\n" +
					"");
		}else{
			int binW = Args.parseInteger(args,"binwidth", 50);
			int binS = Args.parseInteger(args,"binstep", 25);
			ExperimentManager manager = new ExperimentManager(config);
			
			BackgroundDetector detector = new BackgroundDetector(config, manager, binW, binS);
			detector.execute();
			detector.print();
			
			manager.close();
		}
	}

}
