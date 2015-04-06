package edu.psu.compbio.seqcode.projects.seed;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.deepseq.stats.BackgroundCollection;
import edu.psu.compbio.seqcode.deepseq.stats.PoissonBackgroundModel;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.projects.seed.features.EnrichedFeature;
import edu.psu.compbio.seqcode.projects.seed.features.Feature;


/**
 * FeatureDetection is the parent class for all enrichment-based peak/domain callers in this package. 
 *  
 * Handles genome & experiment & configuration loading. 
 * 
 * @author shaunmahony
 *
 */
public abstract class FeatureDetection {
	
	protected GenomeConfig gconfig;
	protected Genome gen;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	protected SEEDConfig sconfig;
	protected Map<Sample, BackgroundCollection> sampleBackgrounds;
	protected Map<ExperimentCondition, BackgroundCollection> conditionBackgrounds;
	protected Map<ExperimentCondition, List<Feature>> features; //The discovered features (defined per condition)
	
	protected boolean strandedEventDetection = false; //Strand aware event detection?
	
	
	/**
	 * Constructor
	 * @param args : command-line arguments
	 */
	public FeatureDetection(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, ExperimentManager man){
		gconfig = gcon;
		econfig = econ;
		sconfig = scon;
		manager = man;
		gen = gconfig.getGenome();
		
		EnrichedFeature.manager = manager;
		
		features=new HashMap<ExperimentCondition, List<Feature>>();
		for(ExperimentCondition c : manager.getConditions())
    		features.put(c,  new ArrayList<Feature>());
		
		sconfig.makeSEEDOutputDirs();
		
		initializeBackgrounds();
	}
	
	/**
	 * Return the full class name of the SEED implementation
	 * @return
	 */
	public abstract String getProgramName();
	
	/**
	 * Return an instance of the implementation-specific thread for this region
	 * @param r : Region
	 */
	public abstract FeatureDetectionThread getMyThread(List<Region> regs);
	
	/**
	 * Post-process the features, and print output. 
	 * Obviously, assume that all features have been collected. 
	 * @return
	 */
	public abstract Map<ExperimentCondition, List<Feature>> postProcess();
	
	/**
	 * The execute method should instantiate the FeatureDetection threads to perform per-region analysis,
	 * and should then collect the Features discovered in each thread for any necessary post-processing. 
	 * This method should also take care of any output printing, etc. 
	 * 
	 * @return : Lists of final Features
	 */
	public Map<ExperimentCondition, List<Feature>> execute(){
		//Split the jobs into the allowed number of threads
		Iterator<Region> testRegions = new ChromosomeGenerator().execute(gen);
		//Threading divides analysis over entire chromosomes. This approach is not compatible with file caching. 
		int numThreads = econfig.getCacheAllData() ? sconfig.getMaxThreads() : 1;
				
		Thread[] threads = new Thread[numThreads];
        List<Region> threadRegions[] = new ArrayList[numThreads];
        int i = 0;
        for (i = 0 ; i < threads.length; i++) {
            threadRegions[i] = new ArrayList<Region>();
        }
        while(testRegions.hasNext()){
        	Region r = testRegions.next(); 
            threadRegions[(i++) % numThreads].add(r);
        }

        for (i = 0 ; i < threads.length; i++) {
        	//Implementation-specific part is in thread's findFeatures method
            Thread t = new Thread(getMyThread(threadRegions[i]));
            t.start();
            threads[i] = t;
        }
        boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(2000);
            } catch (InterruptedException e) { }
            for (i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }
        }
        //Sort the features
        for(ExperimentCondition c : manager.getConditions())
        	Collections.sort(features.get(c));
		
        //Implementation-specific part is in postProcess
        postProcess();
        
        return features;
	}
	
	
	/* Java doesn't allow me to impose this, but it would also be good practise to include in each subclass a 
	 * static method that returns a String describing the command-line arguments. 
	 * The main method in the subclass can then draw on that String if no args are provided. 
	 */
	
	/**
	 * Accessor for the list of features discovered for a given condition.
	 * @param c
	 * @return
	 */
	public List<Feature> getFeatures(ExperimentCondition c){
		if(features==null)
			return null;
		else if(!features.containsKey(c))
			return new ArrayList<Feature>();
		else
			return features.get(c);
	}
	
	/**
	 * Print the features to output files. 
	 * Prefix is defined by SEEDConfig. Suffix is defined depending on what the lists are. 
	 * @param f
	 */
	protected void printEventsFile(Map<ExperimentCondition, List<Feature>> fLists, String suffix){
		try {
			for(ExperimentCondition cond : manager.getConditions()){
				String condName = cond.getName().equals("experiment") ? "" : "_"+cond.getName();
	    		String filename = sconfig.getOutputParentDir()+File.separator+sconfig.getOutBase()+condName+suffix;
				FileWriter fout = new FileWriter(filename);
				//Headers
				fout.write("#"+getProgramName()+"\n");
				fout.write("#Arguments:\t"+sconfig.getArgs()+"\n");
				fout.write("##date "+(new Date()).toString()+"\n");
				//Features
				int fcount=0;
				for(Feature f : fLists.get(cond)){
					if(fcount==0 && !sconfig.outputGFF())
						fout.write(f.headString()+"\n");
					
					if(sconfig.outputGFF())
						fout.write(f.toGFF()+"\n");
					else
						fout.write(f.toString()+"\n");
					fcount++;
		    	}
				fout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Set up the background models
	 * 
	 * By default, the background collections contain only a global Poisson model. 
	 * If the use of local backgrounds is requested, the following occurs for each CONDITION backgorund model:
	 * 		If there is a control and the control is scaled, a local threshold based on expectation from CONTROL is added to the signal
	 * 		If there is no control or the control is not scaled, a local threshold based on expectation from SIGNAL is added to the signal
	 * This makes the behavior of local background thresholds similar to that used by MACS. 
	 */
	protected void initializeBackgrounds(){
		sampleBackgrounds = new HashMap<Sample, BackgroundCollection>();
		conditionBackgrounds = new HashMap<ExperimentCondition, BackgroundCollection>();
		
		//Set up data structures
		for(Sample s : manager.getSamples())
			sampleBackgrounds.put(s, new BackgroundCollection());
		for(ExperimentCondition c : manager.getConditions())
			conditionBackgrounds.put(c, new BackgroundCollection());
		
        //Non-stranded
        if(!strandedEventDetection){
    		//Sample-level genomic backgrounds
        	for(Sample s : manager.getSamples())
        		sampleBackgrounds.get(s).addBackgroundModel(new PoissonBackgroundModel(-1, sconfig.getPerBinPoissonLogPThres(), s.getHitCount(), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '.', 1, true));
        	//Condition-level genomic backgrounds
        	for(ExperimentCondition c : manager.getConditions())
        		conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(-1, sconfig.getPerBinPoissonLogPThres(), c.getTotalSignalCount(), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '.', 1, true));
        	
            //Condition-level locals backgrounds
            for(Integer i : sconfig.getLocalBackWins()){
            	for(ExperimentCondition c : manager.getConditions()){
            		if(c.getControlSamples().size()>0){
            			conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), sconfig.getPerBinPoissonLogPThres(), c.getTotalControlCount(), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '.', c.getPooledSampleControlScaling(), false));
            		}else{
            			//local signal high -- this may bias against locally enriched signal regions, and so should only be used if there is no control or if the control is not yet scaled
                    	if(i.intValue()>=5000) // we don't want the window too small in this case
                    		conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), sconfig.getPerBinPoissonLogPThres(), c.getTotalSignalCount(), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '.', c.getPooledSampleControlScaling(), true));
            		}
            	}
            }            
        }else{ //Stranded
        	//Sample-level genomic backgrounds
        	for(Sample s : manager.getSamples()){
        		sampleBackgrounds.get(s).addBackgroundModel(new PoissonBackgroundModel(-1, sconfig.getPerBinPoissonLogPThres(), s.getStrandedHitCount('+'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '+', 1, true));
        		sampleBackgrounds.get(s).addBackgroundModel(new PoissonBackgroundModel(-1, sconfig.getPerBinPoissonLogPThres(), s.getStrandedHitCount('-'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '-', 1, true));
        	}
        	//Condition-level genomic backgrounds
        	for(ExperimentCondition c : manager.getConditions()){
        		conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(-1, sconfig.getPerBinPoissonLogPThres(), c.getStrandedTotalSignalCount('+'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '+', 1, true));
        		conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(-1, sconfig.getPerBinPoissonLogPThres(), c.getStrandedTotalSignalCount('-'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '-', 1, true));
        	}
        	
            //Condition-level locals backgrounds
            for(Integer i : sconfig.getLocalBackWins()){
            	for(ExperimentCondition c : manager.getConditions()){
            		if(c.getControlSamples().size()>0){
            			conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), sconfig.getPerBinPoissonLogPThres(), c.getStrandedTotalControlCount('+'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '+', c.getPooledSampleControlScaling(), false));
            			conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), sconfig.getPerBinPoissonLogPThres(), c.getStrandedTotalControlCount('-'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '-', c.getPooledSampleControlScaling(), false));
            		}else{
            			//local signal high -- this may bias against locally enriched signal regions, and so should only be used if there is no control or if the control is not yet scaled
                    	if(i.intValue()>=5000){ // we don't want the window too small in this case
                    		conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), sconfig.getPerBinPoissonLogPThres(), c.getStrandedTotalSignalCount('+'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '+', c.getPooledSampleControlScaling(), true));
                    		conditionBackgrounds.get(c).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), sconfig.getPerBinPoissonLogPThres(), c.getStrandedTotalSignalCount('-'), gen.getGenomeLength(), econfig.getMappableGenomeProp(), (double)sconfig.getBinWidth(), '-', c.getPooledSampleControlScaling(), true));
                    	}
            		}
            	}
            }
        }
	}
	
	/**
	 * Filter the features by score & size
	 * @param features
	 * @param threshold
	 * @param below : accept only features below threshold
	 * @return
	 */
	protected Map<ExperimentCondition, List<Feature>> filter(Map<ExperimentCondition, List<Feature>> features, double threshold, boolean below){
		Map<ExperimentCondition, List<Feature>> results = new HashMap<ExperimentCondition, List<Feature>>();
		for(ExperimentCondition cond : manager.getConditions()){
			List<Feature> cres = new ArrayList<Feature>(); 
			for(Feature f : features.get(cond)){
				if(f.getCoords().getWidth()>= sconfig.getMinFeatureSize()){
					if(below && f.getScore()<threshold)
						cres.add(f);
					else if(!below && f.getScore()>threshold)
						cres.add(f);
				}
			}
			results.put(cond, cres);
		}
		return results;
	}
	
	/**
	 * FeatureDetectionThread: thread that searches for the features of interest
	 * @author mahony
	 *
	 */
	public abstract class FeatureDetectionThread implements Runnable {
        protected Collection<Region> runRegions;
        //Hits maintained in separate lists per strands - it's easier to do feature trimming & quantification this way 
        protected Map<Sample, List<StrandedBaseCount>> hitsPos; 	//Lists of positive strand tags in the current region. Indexed by Sample.
        protected Map<Sample, List<StrandedBaseCount>> hitsNeg; 	//Lists of negative strand tags in the current region. Indexed by Sample.
        protected float[][][] landscape=null;  		//Binned tag density in the current region after shifting and extending. Indexed by Sample, base, strand
        protected Map<ExperimentCondition, List<Feature>> threadFeatures;
        protected int shift=sconfig.getTagShift(), hit3Extend=sconfig.getTag3PrimeExtension(), hit5Extend=sconfig.getTag5PrimeExtension();
        
        public FeatureDetectionThread(List<Region> regs){
        	runRegions = regs;
        	threadFeatures = new HashMap<ExperimentCondition, List<Feature>>();
        	for(ExperimentCondition c : manager.getConditions())
        		threadFeatures.put(c,  new ArrayList<Feature>());
        }
        /**
         * Run the thread, executing feature detection on each listed region.
         * @param regs : regions that this thread examines
         */
		public void run(){ //declare as final if not overloaded (guarantees that landscapes is correctly set)
			//Run over each region
			for (Region currentRegion : runRegions) {
                //Split the job up into large chunks
                for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=sconfig.MAXSECTION){
                    int y = (int) (x+sconfig.MAXSECTION);
                    if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                    Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
                    
                    hitsPos = new HashMap<Sample, List<StrandedBaseCount>>();
                    hitsNeg = new HashMap<Sample, List<StrandedBaseCount>>();
                
                    //Initialize & sort the read lists per Sample
                	for(Sample samp : manager.getSamples()){
                		synchronized(manager){//hitCache requires thread safety
                			List<StrandedBaseCount> sampHitsP = samp.getStrandedBases(currSubRegion, '+'); 
                			List<StrandedBaseCount> sampHitsN = samp.getStrandedBases(currSubRegion, '-');
                			Collections.sort(sampHitsP); Collections.sort(sampHitsN); //This might be pointless - the hits should be sorted in the cache already
                			hitsPos.put(samp, sampHitsP);
                			hitsNeg.put(samp, sampHitsN);
                		}
                	}
                	//makeHitLandscape & make GaussianLandscape populate the landscape data structure
                	//execute can therefore assume that these structures are updated, unless run() is overloaded
                    if(sconfig.getBinWidth()==1 && sconfig.getTagGaussSigma()>0)
                    	makeGaussianLandscape(hitsPos, hitsNeg, currSubRegion, sconfig.getTagGaussSigma(), sconfig.getTagGaussWidth());
                    else
                    	makeHitLandscape(hitsPos, hitsNeg, currSubRegion, sconfig.getBinWidth(), sconfig.getBinStep());
                    
                    //Implementation-specific execution
                    Map<ExperimentCondition, List<Feature>> currFeatures = findFeatures(currSubRegion);
                    
                    //Add to thread's features
                    for(ExperimentCondition cond : manager.getConditions())
                    	threadFeatures.get(cond).addAll(currFeatures.get(cond));
                    
                    //TODO: Progress print (replace with some kind of percentage update)
        			System.err.print(".");
                }
			}
			//Filter excluded regions (if necessary)
			threadFeatures = filterExcluded(threadFeatures);
			
			//Add all threadFeatures to the overall results
			synchronized(features){
				for(ExperimentCondition cond : manager.getConditions())
					features.get(cond).addAll(threadFeatures.get(cond));
			}
		}
		
		/**
		 * The core functionality in any event finder should be implemented in this method.
		 * Assumes hits, landscape has been initialized
		 * @param r : Region to analyze
		 * @return : Lists of Features in each ExperimentCondition in the current region
		 */
		public abstract Map<ExperimentCondition, List<Feature>> findFeatures(Region r);
		
		
		/**
         * Makes integer arrays corresponding to the read landscape over the current region.
         * Tags are semi-extended by half the binWidth to account for the step, 
         * 	  and may also be shifted or extended here, depending on the event detection strategy
         * No needlefiltering here as that is taken care of during tag loading (i.e. in Sample)
         * 
         * @param hits  : Lists of StrandedBaseCounts, indexed by Sample (sorted within each sample)
         * @param currReg
         * @param binWidth
         * @param binStep
         */
    	protected void makeHitLandscape(Map<Sample, List<StrandedBaseCount>> hitsPos, Map<Sample, List<StrandedBaseCount>> hitsNeg, Region currReg, int binWidth, int binStep){
    		int numBins = (int)(currReg.getWidth()/binStep);
    		landscape = new float[hitsPos.size()][numBins+1][2];
    		int halfWidth = binWidth/2;

    		for(Sample samp : manager.getSamples()){
    			for(int i=0; i<=numBins; i++)
    				for(int s=0; s<=1; s++)
    					landscape[samp.getIndex()][i][s]=0;
    			
    			for(int strand=0; strand<=1; strand++){
	        		List<StrandedBaseCount> currHits = strand==0 ? hitsPos.get(samp) : hitsNeg.get(samp);
	    			
		    		for(StrandedBaseCount h : currHits){
		    			
		    			//landscape array
		    			int left = getLeft(h);
		    			int right = getRight(h);
		    			if(left <= currReg.getEnd() && right>=currReg.getStart()){
			    			int offsetL=inBounds(left-currReg.getStart(),0,currReg.getWidth());
			    			int offsetR=inBounds(right-currReg.getStart(),0,currReg.getWidth());
		    				
			    			int binstart = inBounds(((offsetL-halfWidth)/binStep), 0, numBins);
		    				int binend = inBounds(((offsetR/binStep)), 0, numBins);
		    				for(int b=binstart; b<=binend; b++)
		    					landscape[samp.getIndex()][b][strand]+=h.getCount();
		    			}
	            	}
    			}
    		}
    	}

    	/**
    	 * Makes integer arrays corresponding to the Gaussian-smoothed read landscape over the current region.
    	 * This only operates at single-bp resolution.
         * Tags are gaussian smoothed over the landscape, and may also be shifted, depending on the event detection strategy
         * No needlefiltering here as that is taken care of during tag loading (i.e. in Sample)
         * 
         * @param hits  : Lists of StrandedBaseCounts, indexed by Sample
         * @param currReg
         * @param gaussSigma: Gaussian sigma (std dev)
         * @param gaussWidth: width over which to 'extend' each tag
    	 */
    	protected void makeGaussianLandscape(Map<Sample, List<StrandedBaseCount>> hitsPos, Map<Sample, List<StrandedBaseCount>> hitsNeg, Region currReg, float gaussSigma, int gaussWidth){
    		int length = (int)currReg.getWidth();
    		landscape = new float[hitsPos.size()][length+1][2];
    		float [][][] fivePrimes = new float[hitsPos.size()][length+1][2];
    		float[] kernel = initGaussianKernel(gaussSigma, gaussWidth);
    		
    		for(Sample samp : manager.getSamples()){
    			for(int i=0; i<=length; i++)
    				for(int s=0; s<=1; s++){
    					landscape[samp.getIndex()][i][s]=0; fivePrimes[samp.getIndex()][i][s]=0;
    				}
    			for(int strand=0; strand<=1; strand++){
	        		List<StrandedBaseCount> currHits = strand==0 ? hitsPos.get(samp) : hitsNeg.get(samp);
	    			
		    		for(StrandedBaseCount h : currHits){
		    			//(shifted) fivePrimes array
		    			int offset5=inBounds(getShifted5Prime(h)-currReg.getStart(),0,currReg.getWidth());
		    			int binoff5 = inBounds((int)(offset5), 0, length);
		    			fivePrimes[samp.getIndex()][binoff5][strand]+=h.getCount();
		    		}
		    		//landscape array is fivePrime * gaussian 
		    		float total=0;
		    		for(int s=0; s<=1; s++){
			    		for (int i=0;i<length;i++){
			    			float v=kernel[0]*fivePrimes[samp.getIndex()][i][s] + Float.MIN_VALUE;		// init with very small number
			                float weight=kernel[0];
			                for (int j = 1; j < kernel.length && i+j < length; j++) {
			                    v+=fivePrimes[samp.getIndex()][i+j][s]*kernel[j];
			                    weight += kernel[j];                
			                }
			                for (int j = 1; j < kernel.length && i-j >= 0; j++) {
			                    v+=fivePrimes[samp.getIndex()][i-j][s]*kernel[j];
			                    weight += kernel[j];                
			                }
			    			v = v / weight;
			    			landscape[samp.getIndex()][i][s] = v;
			    			total+=v;
			    		}
			    		for (int i=0;i<length;i++)
			    			landscape[samp.getIndex()][i][s]=landscape[samp.getIndex()][i][s]/total;
		    		}
    			}
    		}
    	}
    	
    	protected final int inBounds(int x, int min, int max){
    		if(x<min){return min;}
    		if(x>max){return max;}
    		return x;
    	}
    	protected final int getShifted5Prime(StrandedBaseCount h){
    		return(h.getStrand()=='+' ? 
    				h.getCoordinate()+shift : 
    				h.getCoordinate()-shift);
    	}
    	protected final int getLeft(StrandedBaseCount h){
    		return(h.getStrand()=='+' ? 
    				h.getCoordinate()+shift-hit5Extend : 
    				h.getCoordinate()-shift-hit3Extend);
    	}
    	protected final int getRight(StrandedBaseCount h){
    		return(h.getStrand()=='+' ? 
    				h.getCoordinate()+shift+hit3Extend : 
    				h.getCoordinate()-shift+hit5Extend);
    	}
    	/**
    	 * Initializes half of a Gaussian
    	 * @param gaussSigma
    	 * @param gaussWidth
    	 * @return
    	 */
    	protected float[] initGaussianKernel(float gaussSigma, int gaussWidth){
    		float[] y = new float[(gaussWidth/2)+1];
    		//gaussWidth*gaussWidth = variance
    		NormalDistribution gaussian = new NormalDistribution(0, gaussWidth*gaussWidth);
    		float total=0;
    		for(int i=0; i<y.length; i++){
    			y[i] = (float)gaussian.calcProbability((double)i);
    			total += y[i];
    		}
    		//Normalize
    		for(int i=0; i<y.length; i++)
    			y[i] /= total;
    		
    		return y;
    	}
    	
    	/**
    	 * Parses the landscape arrays to get a per-condition count array 
    	 * 
    	 * @param cond : ExperimentCondition of interest
    	 * @param data : data structure to parse. Should be landscape. Assumes indexed by Sample, base, strand
    	 * @param strand : +/-/.
    	 * @param signal : true to count condition's signal samples, false to count condition's control samples. 
    	 * @return
    	 */
    	protected float[]  getConditionCounts(ExperimentCondition cond, float[][][] data, char strand, boolean signal){
    		int clength = data[0].length;
    		float[] counts = new float[clength];
    		for(int c=0; c<clength; c++){
	    		float currCount=0;
	    		List<Sample> currSamples = signal ? cond.getSignalSamples() : cond.getControlSamples();
				for(Sample samp : currSamples){
					if(strand=='.' || strand=='+')
						currCount+=data[samp.getIndex()][c][0];
					if(strand=='.' || strand=='-')
						currCount+=data[samp.getIndex()][c][1];
				}
				counts[c]=currCount;
    		}
			return counts;
    	}
    	
    	/**
    	 * Uses a rough, inexact binary search to find StrandedBaseCounts that overlap a given feature in each Sample
    	 * @param sameStrHits : Lists of StrandedBaseCounts, all from same strand, indexed by Sample - assumes sorted
    	 * @param f : Feature
    	 * @return : Lists of StrandedBaseCounts that overlap the feature coordinates
    	 */
    	protected Map<Sample, List<StrandedBaseCount>> overlappingHits(Map<Sample, List<StrandedBaseCount>> sameStrHits, Feature f){
    		Map<Sample, List<StrandedBaseCount>> subHits = new HashMap<Sample, List<StrandedBaseCount>>();
    		for(Sample samp : manager.getSamples()){
    			List<StrandedBaseCount> sub = new ArrayList<StrandedBaseCount>();
	    		int l = 0;
	            int r = sameStrHits.get(samp).size();
	            int featureStart = f.getCoords().getStart();
	            int featureEnd= f.getCoords().getEnd();
	            char str = f.getCoords().getStrand();
	            while (r - l > 10) {
	                int c = (l + r) / 2;
	                if (featureStart > getRight(sameStrHits.get(samp).get(c))) {
	                    l = c;
	                } else {
	                    r = c;
	                }
	            }
	            while (l > 0 && (getRight(sameStrHits.get(samp).get(l)) >= featureStart)) {
	                l--;
	            }
	            while (l < sameStrHits.get(samp).size() && getLeft(sameStrHits.get(samp).get(l)) <= featureEnd){
	            	StrandedBaseCount hit = sameStrHits.get(samp).get(l);
	            	int hitL = getLeft(hit);
	            	int hitR = getRight(hit);
	    			if(f.getCoords().overlaps(hitL, hitR) 
	    					&& (str=='.' || str==hit.getStrand())){
	    				sub.add(hit);
	    			}			
	                l++;
	            }
	            subHits.put(samp, sub);
    		}
    		return(subHits);
    	}
    	
    	/**
    	 * Filter out pre-defined regions to ignore (e.g. tower / blacklist regions)
    	 * @param testRegions
    	 * @return
    	 */
        protected Map<ExperimentCondition, List<Feature>> filterExcluded(Map<ExperimentCondition, List<Feature>> testFeatures) {
        	if(sconfig.getRegionsToIgnore().size()==0)
    			return testFeatures;
        	
        	for (ExperimentCondition ec : manager.getConditions()){
				for(Region blackList : sconfig.getRegionsToIgnore()){
					Iterator<Feature> it = testFeatures.get(ec).iterator();
					while(it.hasNext()){
						Feature ef = it.next();
						if(ef.getCoords().overlaps(blackList)){
							it.remove();
						}
					}
				}
			}
        	return testFeatures;
    	}
    	
	}
}
