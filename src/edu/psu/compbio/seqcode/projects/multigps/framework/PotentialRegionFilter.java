package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;

/**
 * PotentialRegionFilter: Find a set of regions that are above a threshold in at least one replicate. 
 * 		A region the size of the model span (i.e. 2x model range) potentially contains a binding site if 
 * 		it passes all Poisson thresholds in at least one replicate from one condition.
 * 		The Poisson thresholds are based on the model span size to keep consistent with the final used thresholds. 
 * Overall counts for reads in potential regions and outside potential regions are maintained to assist noise model initialization.  
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class PotentialRegionFilter {

	protected ExperimentManager manager; 
	protected Config config;
	protected Genome gen;
	protected float maxBinWidth=0, binStep, winExt;
	protected boolean loadControl=true; 
	protected boolean stranded=false;
	protected List<Region> potentialRegions = new ArrayList<Region>();
	protected double potRegionLengthTotal=0;
	protected HashMap<ControlledExperiment, BackgroundCollection> replicateBackgrounds=new HashMap<ControlledExperiment, BackgroundCollection>(); //Background models for each replicate
	protected HashMap<ControlledExperiment, Double> potRegCountsSigChannel = new HashMap<ControlledExperiment, Double>();
	protected HashMap<ControlledExperiment, Double> nonPotRegCountsSigChannel = new HashMap<ControlledExperiment, Double>();
	protected HashMap<ControlledExperiment, Double> potRegCountsCtrlChannel = new HashMap<ControlledExperiment, Double>();
	protected HashMap<ControlledExperiment, Double> nonPotRegCountsCtrlChannel = new HashMap<ControlledExperiment, Double>();
	
	public PotentialRegionFilter(Config c, ExperimentManager man){
		manager = man;
		config = c; 
		gen = config.getGenome();
		//Initialize background models
		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
    		for(ControlledExperiment rep : cond.getReplicates()){
    			float binWidth = rep.getBindingModel().getInfluenceRange();
    			if(binWidth>maxBinWidth){maxBinWidth=binWidth;}
    			replicateBackgrounds.put(rep, new BackgroundCollection());
    			//global threshold
    			replicateBackgrounds.get(rep).addBackgroundModel(new PoissonBackgroundModel(-1, config.getPRLogConf(), (1-rep.getSigProp())*rep.getSignal().getHitCount(), config.getGenome().getGenomeLength(), config.getMappableGenomeProp(), rep.getBindingModel().getInfluenceRange(), '.', 1, true));
    			for(Integer i : config.getLocalBackgroundWindows()){
                	if(rep.hasControl()){//local control thresholds 
                    	//signal threshold based on what would be expected from the CONTROL locality
                		replicateBackgrounds.get(rep).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), config.getPRLogConf(), rep.getControl().getHitCount(),  config.getGenome().getGenomeLength(), config.getMappableGenomeProp(), rep.getBindingModel().getInfluenceRange(), '.', rep.getControlScaling(), false));
                	}else{
                		//local signal threshold -- this may bias against locally enriched signal regions, and so should only be used if there is no control or if the control is not yet scaled
                    	if(i.intValue()>=10000) // we don't want the window too small in this case
                    		replicateBackgrounds.get(rep).addBackgroundModel(new PoissonBackgroundModel(i.intValue(), config.getPRLogConf(), (1-rep.getSigProp())*rep.getSignal().getHitCount(), config.getGenome().getGenomeLength(), config.getMappableGenomeProp(), rep.getBindingModel().getInfluenceRange(), '.', 1, true));
                	}	
                }
    			
    			System.err.println("PotentialRegionFilter: genomic threshold for "+rep.getName()+" with bin width "+rep.getBindingModel().getInfluenceRange()+" = "+replicateBackgrounds.get(rep).getGenomicModelThreshold());
    			
    			//Initialize counts
    			potRegCountsSigChannel.put(rep, 0.0);
    			nonPotRegCountsSigChannel.put(rep, 0.0);
    			potRegCountsCtrlChannel.put(rep, 0.0);
    			nonPotRegCountsCtrlChannel.put(rep, 0.0);
    		}
    	}
		binStep = config.POTREG_BIN_STEP;
		if(binStep>maxBinWidth/2)
			binStep=maxBinWidth/2;
		winExt = maxBinWidth/2;
	}
	
	//Accessors for read counts
	public Double getPotRegCountsSigChannel(ControlledExperiment e){ return potRegCountsSigChannel.get(e);}
	public Double getNonPotRegCountsSigChannel(ControlledExperiment e){ return nonPotRegCountsSigChannel.get(e);}
	public Double getPotRegCountsCtrlChannel(ControlledExperiment e){ return potRegCountsCtrlChannel.get(e);}
	public Double getNonPotRegCountsCtrlChannel(ControlledExperiment e){ return nonPotRegCountsCtrlChannel.get(e);}
	public List<Region> getPotentialRegions(){return potentialRegions;}
	public double getPotRegionLengthTotal(){return potRegionLengthTotal;}
	
	/**
	 * Find list of potentially enriched regions 
	 * (windows that contain the minimum number of reads needed to pass the Poisson backgrounds).
	 * @param testRegions
	 */
	public List<Region> execute(){
		//TODO: check config for defined subset of regions
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
        
		//Initialize signal & noise counts based on potential region calls
        for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
    		for(ControlledExperiment rep : cond.getReplicates()){
    			if(rep.getSigProp()==0) //Only update if not already initialized
    				rep.setSigNoiseCounts(potRegCountsSigChannel.get(rep), nonPotRegCountsSigChannel.get(rep));
    		}
        }
        
        for(Region r : potentialRegions)
        	potRegionLengthTotal+=(double)r.getWidth();
        
     	return potentialRegions;
	}
	
	/**
     * Print potential regions to a file.
     * TESTING ONLY 
     */
    public void printPotentialRegionsToFile(){
    	try {
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".potential.regions";
			FileWriter fout = new FileWriter(filename);
			for(Region r : potentialRegions){
	    		fout.write(r.getLocationString()+"\n");			
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
	
    class PotentialRegionFinderThread implements Runnable {
        private Collection<Region> regions;
        private double[][] landscape=null;
        private double[][] starts=null;
        private List<Region> threadPotentials = new ArrayList<Region>();
        
        public PotentialRegionFinderThread(Collection<Region> r) {
            regions = r;
        }
        
        public void run() {
        	int expansion = (int)(winExt + maxBinWidth/2);
        	for (Region currentRegion : regions) {
            	Region lastPotential=null;
                //Split the job up into large chunks
                for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
                    int y = (int) (x+config.MAXSECTION+(expansion)); //Leave a little overhang to handle enriched regions that may hit the border. Since lastPotential is defined above, a region on the boundary should get merged in.
                    if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                    Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
                    
                    List<Region> currPotRegions = new ArrayList<Region>();
                    List<List<StrandedBaseCount>> ipHits = new ArrayList<List<StrandedBaseCount>>();
                    List<List<StrandedBaseCount>> backHits = new ArrayList<List<StrandedBaseCount>>();
                    
                    synchronized(manager){
	                    //Initialize the read lists
                    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
                    		for(ControlledExperiment rep : cond.getReplicates()){
                    			ipHits.add(new ArrayList<StrandedBaseCount>());
                    			backHits.add(new ArrayList<StrandedBaseCount>());
                    		}
                    	}
                    	//Load reads by replicate
                    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
                    		for(ControlledExperiment rep : cond.getReplicates()){
                    			ipHits.get(rep.getIndex()).addAll(rep.getSignal().getUnstrandedBases(currSubRegion));
                    			if(loadControl && rep.hasControl())
                    				backHits.get(rep.getIndex()).addAll(rep.getControl().getUnstrandedBases(currSubRegion));
                    		}
                    	}
                    }
            		int numStrandIter = stranded ? 2 : 1;
                    for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
                        //If stranded peak-finding, run over both strands separately
                        char str = !stranded ? '.' : (stranditer==1 ? '+' : '-');
					 
                        makeHitLandscape(ipHits, currSubRegion, maxBinWidth, binStep, str);
                        double ipHitCounts[][] = landscape.clone();
                        double ipBinnedStarts[][] = starts.clone();
                        //double backHitCounts[][] = null;
                        double backBinnedStarts[][] = null;
                        if (loadControl) {
                            makeHitLandscape(backHits, currSubRegion, maxBinWidth, binStep, str);
                            //backHitCounts = landscape.clone();
                            backBinnedStarts = starts.clone();
                        }
					
                        //Scan regions
                        int currBin=0;
                        for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)maxBinWidth; i+=(int)binStep){
                        	boolean regionPasses=false;
                        	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
                            	for(ControlledExperiment rep : cond.getReplicates()){
		                            double ipWinHits=ipHitCounts[rep.getIndex()][currBin];
		                            //First Test: is the read count above the genome-wide thresholds? 
		                            if(replicateBackgrounds.get(rep).passesGenomicThreshold((int)ipWinHits, str)){
		                            	//Second Test: refresh all thresholds & test again
		                            	replicateBackgrounds.get(rep).updateModels(currSubRegion, i-x, ipBinnedStarts[rep.getIndex()], backBinnedStarts==null ? null : backBinnedStarts[rep.getIndex()], binStep);
		                            	if(replicateBackgrounds.get(rep).passesAllThresholds((int)ipWinHits, str)){
		                                    //If the region passes the thresholds for one replicate, it's a potential
		                            		regionPasses=true;
		                            		break;
		                                }
		                            }
		                        }
                        	}
                        	if(regionPasses){
                        		Region currPotential = new Region(gen, currentRegion.getChrom(), Math.max(i-expansion, 1), Math.min((int)(i-1+expansion), currentRegion.getEnd()));
                        		if(lastPotential!=null && currPotential.overlaps(lastPotential)){
                        			lastPotential = lastPotential.expand(0, currPotential.getEnd()-lastPotential.getEnd());
                        		}else{
                        			//Add the last recorded region to the list
                        			if(lastPotential!=null){
                        				if(lastPotential.getWidth()<=config.getBMAnalysisWindowMax()){
                        					currPotRegions.add(lastPotential);
                        					threadPotentials.add(lastPotential);
                        				}else{
                        					//Break up long windows
                        					List<Region> parts = breakWindow(lastPotential, ipHits, config.getBMAnalysisWindowMax(), str);
                        					for(Region p : parts){
                        						currPotRegions.add(p);
                        						threadPotentials.add(p);
                        					}
                        				}
                        			}lastPotential = currPotential;
                        		}
                        	}
                            currBin++;
                        }
					}
                    //Count all "signal" reads overlapping the regions in currPotRegions (including the lastPotential)
                    if(lastPotential!=null)
                    	currPotRegions.add(lastPotential);
                    currPotRegions = filterExcluded(currPotRegions);
                    countReadsInRegions(currPotRegions, ipHits, backHits, y==currentRegion.getEnd() ? y : y-expansion);
                    //Note: it looks like currPotRegions and threadPotentials are redundant in the above, but they are not.
                    //currPotRegions is only used to count sig/noise reads in the current section. threadPotentials stores regions over the entire run.
                }
                //Add the final recorded region to the list
                if(lastPotential!=null)
    				threadPotentials.add(lastPotential);
                threadPotentials = filterExcluded(threadPotentials);
            }
        	if(threadPotentials.size()>0){
        		synchronized(potentialRegions){
        			potentialRegions.addAll(threadPotentials);
        		}
        	}	
        }
        
        //Break up a long window into parts
        //For now, we just choose the break points as the bins with the lowest total signal read count around the desired length.
        //TODO: improve?
        protected List<Region> breakWindow(Region lastPotential, List<List<StrandedBaseCount>> ipHits, int preferredWinLen, char str) {
			List<Region> parts = new ArrayList<Region>();
			makeHitLandscape(ipHits, lastPotential, maxBinWidth, binStep, str);
            double ipHitCounts[][] = landscape.clone();
            
            int currPartStart = lastPotential.getStart();
            double currPartTotalMin=Double.MAX_VALUE; int currPartTotalMinPos = -1;
            int currBin=0;
            for(int i=lastPotential.getStart(); i<lastPotential.getEnd()-(int)maxBinWidth; i+=(int)binStep){
            	if(lastPotential.getEnd()-currPartStart < (preferredWinLen*1.5))
            		break;
            	int currBinTotal=0;
            	for(ExperimentCondition cond : manager.getExperimentSet().getConditions())
                	for(ControlledExperiment rep : cond.getReplicates())
                		currBinTotal+=ipHitCounts[rep.getIndex()][currBin];
            	
            	if(i>(currPartStart+preferredWinLen-1000) && i<(currPartStart+preferredWinLen+1000)){ 
            		if(currBinTotal<currPartTotalMin){
            			currPartTotalMin=currBinTotal;
            			currPartTotalMinPos=i;
            		}
            	}
            	//Add a new part
            	if(i>=(currPartStart+preferredWinLen+1000)){
            		parts.add(new Region(lastPotential.getGenome(), lastPotential.getChrom(), currPartStart, currPartTotalMinPos));
            		currPartStart = currPartTotalMinPos+1;
            		currPartTotalMin=Double.MAX_VALUE; currPartTotalMinPos = -1;
            	}
            	currBin++;
            }
            parts.add(new Region(lastPotential.getGenome(), lastPotential.getChrom(), currPartStart, lastPotential.getEnd()));
            
			return parts;
		}

		//Filter out pre-defined regions to ignore (e.g. tower regions)
        protected List<Region> filterExcluded(List<Region> testRegions) {
			List<Region> filtered = new ArrayList<Region>();
			if(config.getRegionsToIgnore().size()==0)
				return testRegions;
			
			for(Region t : testRegions){
				boolean ignore = false;
				for(Region i : config.getRegionsToIgnore()){
					if(t.overlaps(i)){
						ignore = true; break;
					}
				}
				if(!ignore)
					filtered.add(t);
			}
			return filtered;
		}

		//Makes integer arrays corresponding to the read landscape over the current region.
        //Reads are semi-extended out to bin width
        //No needlefiltering here as that is taken care of during read loading (i.e. in Sample)
    	protected void makeHitLandscape(List<List<StrandedBaseCount>> hits, Region currReg, float binWidth, float binStep, char strand){
    		int numBins = (int)(currReg.getWidth()/binStep);
    		landscape = new double[hits.size()][numBins+1];
    		starts = new double[hits.size()][numBins+1];
    		float halfWidth = binWidth/2;

    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
            	for(ControlledExperiment rep : cond.getReplicates()){
            		List<StrandedBaseCount> currHits = hits.get(rep.getIndex());
	    			for(int i=0; i<=numBins; i++){landscape[rep.getIndex()][i]=0; starts[rep.getIndex()][i]=0; }
		    		for(StrandedBaseCount r : currHits){
		    			if(strand=='.' || r.getStrand()==strand){
		    				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
		    				int binoff = inBounds((int)(offset/binStep), 0, numBins);
		    				starts[rep.getIndex()][binoff]+=r.getCount();
		    				int binstart = inBounds((int)((double)(offset-halfWidth)/binStep), 0, numBins);
		    				int binend = inBounds((int)((double)(offset+halfWidth)/binStep), 0, numBins);
		    				for(int b=binstart; b<=binend; b++)
		    					landscape[rep.getIndex()][b]+=r.getCount();
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
    	
    	/**
    	 * Count the total reads within potential regions via semi binary search.
    	 * Assumes both regs and ipHits are sorted.
    	 * We don't have to check chr String matches, as the hits were extracted from the chromosome
    	 * EndCoord accounts for the extra overhang added to some wide regions
    	 * We also ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
    	 * (and that thus will not be accounted for by the global noise model)  
    	 * @param regs
    	 * @param ipHits
    	 * @param ctrlHits
    	 * @param endCoord
    	 */
    	protected void countReadsInRegions(List<Region> regs, List<List<StrandedBaseCount>> ipHits, List<List<StrandedBaseCount>> ctrlHits, int endCoord){
    		//Iterate through experiments
    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
        		for(ControlledExperiment rep : cond.getReplicates()){
        			double currPotWeightSig=0, currNonPotWeightSig=0, currPotWeightCtrl=0, currNonPotWeightCtrl=0;
        			//Iterate through signal hits
        			for(StrandedBaseCount hit : ipHits.get(rep.getIndex())){
        				if(regs.size()==0)
        					currNonPotWeightSig+=hit.getCount();
        				else{
	    					//Binary search for closest region start
	        				int hpoint = hit.getCoordinate();
	        				if(hpoint<endCoord){ //Throw this check in for the overhang
		        				int l = 0, r = regs.size()-1;
		        	            while (r - l > 1) {
		        	                int c = (l + r) / 2;
		        	                if (hpoint >= regs.get(c).getStart()) {
		        	                    l = c;
		        	                } else {
		        	                    r = c;
		        	                }
		        	            }
		        	            boolean inPot = false;
		        	            for(int x=l; x<=r; x++){
		        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
		        	            		currPotWeightSig+=hit.getCount(); inPot=true; break;
		        	            	}
		        	            }
		        	            if(!inPot)
		        	            	currNonPotWeightSig+=hit.getCount();
	        				}
        				}
        			}
        			//Iterate through control hits
        			for(StrandedBaseCount hit : ctrlHits.get(rep.getIndex())){
        				if(regs.size()==0)
        					currNonPotWeightCtrl+=hit.getCount();
        				else{
	        				//Binary search for closest region start
	        				int hpoint = hit.getCoordinate();
	        				if(hpoint<endCoord){ //Throw this check in for the overhang
		        				int l = 0, r = regs.size()-1;
		        	            while (r - l > 1) {
		        	                int c = (l + r) / 2;
		        	                if (hpoint >= regs.get(c).getStart()) {
		        	                    l = c;
		        	                } else {
		        	                    r = c;
		        	                }
		        	            }
		        	            boolean inPot = false;
		        	            for(int x=l; x<=r; x++){
		        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
		        	            		currPotWeightCtrl+=hit.getCount(); inPot=true; break;
		        	            	}
		        	            }
		        	            if(!inPot)
		        	            	currNonPotWeightCtrl+=hit.getCount();
	        				}
        				}
        			}
        			synchronized(potRegCountsSigChannel){
        				potRegCountsSigChannel.put(rep, potRegCountsSigChannel.get(rep)+currPotWeightSig);
        			}
        			synchronized(nonPotRegCountsSigChannel){
        				nonPotRegCountsSigChannel.put(rep, nonPotRegCountsSigChannel.get(rep)+currNonPotWeightSig);
        			}
        			synchronized(potRegCountsCtrlChannel){
        				potRegCountsCtrlChannel.put(rep, potRegCountsCtrlChannel.get(rep)+currPotWeightCtrl);
        			}
        			synchronized(nonPotRegCountsCtrlChannel){
        				nonPotRegCountsCtrlChannel.put(rep, nonPotRegCountsCtrlChannel.get(rep)+currNonPotWeightCtrl);
        			}
        		}
    		}
	    }
    }
    
    
    /**
	 * This main method is only for testing the PotentialRegionFilter
	 * @param args
	 */
	public static void main(String[] args){
		
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("PotentialRegionFilter:");
			System.err.println(config.getArgsList());
		}else{
			ExperimentManager manager = new ExperimentManager(config);
			RealValuedHistogram histo = new RealValuedHistogram(0, 10000, 20);
			
			ExperimentSet eset = manager.getExperimentSet();
			System.err.println("Conditions:\t"+eset.getConditions().size());
			for(ExperimentCondition c : eset.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : eset.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			
			PotentialRegionFilter filter = new PotentialRegionFilter(config, manager);
			List<Region> potentials = filter.execute();
			
			int min = Integer.MAX_VALUE;
			int max = -Integer.MAX_VALUE;
			for(Region r : potentials){
				System.out.println(r.getLocationString()+"\t"+r.getWidth());
				histo.addValue(r.getWidth());
				if(r.getWidth()<min)
					min = r.getWidth();
				if(r.getWidth()>max)
					max = r.getWidth();
			}
			System.out.println("Potential Regions: "+potentials.size());
			System.out.println("Min width: "+min+"\tMax width: "+max);
			histo.printContents();
			//for(ExperimentCondition c : eset.getConditions()){
			//	for(ControlledExperiment r : c.getReplicates()){
			//		System.err.println(c.getName()+"\t"+r.getName()+"\t"+r.getSigCount()+"\t"+r.getNoiseCount());
			//	}
			//}
			
			manager.close();
		}
	}
}
