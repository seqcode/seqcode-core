package edu.psu.compbio.seqcode.projects.multigps.mixturemodel;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.BackgroundCollection;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.PoissonBackgroundModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.PotentialRegionFilter;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;
import edu.psu.compbio.seqcode.projects.multigps.motifs.MotifPlatform;

/**
 * BindingMixture: defines a mixture model over binding events. 
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class BindingMixture {

	protected ExperimentManager manager;
	protected Config config;
	protected PotentialRegionFilter potRegFilter;
	protected List<Region> testRegions;
	protected HashMap<Region, List<List<BindingComponent>>> activeComponents; //Components active after a round of execute()
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds=new HashMap<ExperimentCondition, BackgroundCollection>(); //Genomic Background models for each condition -- used to set alpha values in sparse prior
	protected List<BindingEvent> bindingEvents;
	protected List<Region> regionsToPlot;
	protected int trainingRound=0;
	protected double noisePerBase[];      //Defines global noise 
	protected double relativeCtrlNoise[]; //Defines global noise
	protected HashMap<Region, Double[]> noiseResp = new HashMap<Region, Double[]>(); //noise responsibilities after a round of execute(). Hashed by Region, indexed by condition
	protected MotifPlatform motifFinder;
	
	public BindingMixture(Config c, ExperimentManager eMan, PotentialRegionFilter filter){
		config = c;
		manager = eMan;
		potRegFilter=filter;
		testRegions = filter.getPotentialRegions();
		if(config.getFindingMotifs())
			motifFinder = new MotifPlatform(config, manager, testRegions);
		else 
			motifFinder=null;
		regionsToPlot = config.getRegionsToPlot();
		bindingEvents = new ArrayList<BindingEvent>();
		BindingEvent.setExperimentSet(manager.getExperimentSet());
		BindingEvent.setConfig(config);
		
		activeComponents = new HashMap<Region, List<List<BindingComponent>>>();
		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
			conditionBackgrounds.put(cond, new BackgroundCollection());
			conditionBackgrounds.get(cond).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), cond.getTotalSignalCount(), config.getGenome().getGenomeLength(), config.getMappableGenomeProp(), cond.getMaxModelRange(), '.', 1, true));
			System.err.println("Alpha "+cond.getName()+"\tRange="+cond.getMaxModelRange()+"\t"+(double)conditionBackgrounds.get(cond).getMaxThreshold('.'));
		}
		
		noisePerBase = new double[manager.getNumConditions()];
		relativeCtrlNoise = new double[manager.getNumConditions()];
		initializeGlobalNoise();
		
		
		if(config.useMultiConditionPosPrior()){
			//Calculating the prior variables here for info only -- actual variables calculated in BindingEM. 
			double N = testRegions.size();
			double S = N*config.getProbSharedBinding();
			double L = (double)config.getGenome().getGenomeLength();
			double infoProbAgivenB = Math.log(config.getProbSharedBinding())/Math.log(2);
	        double infoProbAgivenNOTB =  Math.log((N-S)/(L-N+S))/Math.log(2);
			System.err.println("Multi-condition positional priors:\tA given B:"+infoProbAgivenB+"\tA given notB:"+infoProbAgivenNOTB);
		}
	}
	
	
	/**
	 * Initialize the threads and execute the binding mixture
	 * Relies on the following being set: testRegions, bindingModels (from exptMan)
	 * uniformBindingComponents: if true, components are initialized at uniform prob and spacing, 
	 * 					otherwise, components are initialized from activeComponents
	 * EM: if true, run EM, otherwise run ML assignment   
	 */
	public void execute(boolean EM, boolean uniformBindingComponents){
		trainingRound++;
		Thread[] threads = new Thread[config.getMaxThreads()];
        ArrayList<Region> threadRegions[] = new ArrayList[config.getMaxThreads()];
        int i = 0;
        for (i = 0 ; i < threads.length; i++) {
            threadRegions[i] = new ArrayList<Region>();
        }
        for(Region r : testRegions){
            threadRegions[(i++) % config.getMaxThreads()].add(r);
        }

        for (i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new BindingMixtureThread(threadRegions[i], EM, uniformBindingComponents));
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
	}
	
	/**
	 * Return the active components, flattening the hash map first. 
	 * @return
	 */
	public List<List<BindingComponent>> getBindingComponents(){
		List<List<BindingComponent>> comps = new ArrayList<List<BindingComponent>>();
		for(int c=0; c<manager.getNumConditions(); c++)
			comps.add(new ArrayList<BindingComponent>());
		for(Region r : activeComponents.keySet()){
			for(int c=0; c<manager.getNumConditions(); c++){
				comps.get(c).addAll(activeComponents.get(r).get(c));
			}
		}
		return comps;
	}
	
	/**
	 * Return the discovered binding events. Call after ML assignment.
	 * @return
	 */
	public List<BindingEvent> getBindingEvents(){return bindingEvents;}
	
	/**
	 * Return the initialized motif-finder
	 * @return
	 */
	public MotifPlatform getMotifFinder(){return motifFinder;}
	
	/**
     * Update binding models for each replicate given the discovered binding components. 
     * @param left 
     * @param right
     * @param String filename for new distribution file
     * @return double array of log KL values
     */
    public Double[] updateBindingModel(String distribFilename){
    	int left=0,right=0;
    	for(ControlledExperiment rep : manager.getExperimentSet().getReplicates()){
    		BindingModel mod = rep.getBindingModel();
    		int l=-1*mod.getMin(),r=mod.getMax();
    		if(!config.getFixedModelRange()){
    			Pair<Integer, Integer> newEnds = mod.getNewEnds(300,200);
    			l=newEnds.car(); r=newEnds.cdr();
    		}
    		left = Math.min(config.MAX_BINDINGMODEL_WIDTH/2, Math.max(left, l));
    		right = Math.min(config.MAX_BINDINGMODEL_WIDTH/2, Math.max(right, r));
    	}
    	
    	int numReps = manager.getExperimentSet().getReplicates().size();
    	int width = left+right;
    	int readProfileCenter = config.MAX_BINDINGMODEL_WIDTH/2;
    	int offsetLeft = readProfileCenter-left;
		double[][] newModel_plus=new double[numReps][width];
		double[][] newModel_minus=new double[numReps][width];
		double[][] newModel=new double[numReps][width];
		Double[] logKL = new Double[numReps];
		int[] eventCounter = new int[numReps];
		for (int x=0; x<numReps; x++){
			for (int i=0;i<width;i++){
				newModel_plus[x][i]=1;
				newModel_minus[x][i]=1;
			}
			logKL[x]=0.0;
			eventCounter[x]=0;
		}
		
		if(config.doBMUpdate()){
			//Sum read profiles if there are enough binding components
	    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
	    		List<BindingComponent> currComps = new ArrayList<BindingComponent>();
	    		double currAlpha = (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
	    		//Choose which components to include
	    		for(Region r : activeComponents.keySet()){
	    			//1) Include joint events?
	    			if(activeComponents.get(r).get(cond.getIndex()).size()==1 ||(activeComponents.get(r).get(cond.getIndex()).size()>1 && config.getIncludeJointEventsInBMUpdate())){
	    				for(BindingComponent bc : activeComponents.get(r).get(cond.getIndex())){
	    					//2) Component must not be at the edge of the region 
	    					if((bc.getPosition()-r.getStart()>cond.getMaxInfluenceRange()/2) && (r.getEnd()-bc.getPosition()>cond.getMaxInfluenceRange()/2)){
	    						//3) Arbitrary minimum read support for BM components
	    						if(bc.getSumResponsibility()>(config.getMinComponentReadFactorForBM()*currAlpha))
	    							currComps.add(bc);
	    					}
	    				}
	    			}
	    		}
	    		
	    		if (currComps.size()<config.getMinComponentsForBMUpdate()){
	    			System.err.println("The "+cond.getName()+" read distributions cannot be updated due to too few binding components ("+currComps.size()+"<"+config.getMinComponentsForBMUpdate()+").");
	    			for(ControlledExperiment rep : cond.getReplicates())
	    				logKL[rep.getIndex()]=Double.NaN;
	    		}else{
	    			for(ControlledExperiment rep : cond.getReplicates()){
			    		int x = rep.getIndex();
			    		
			    		for(BindingComponent comp : currComps){
				    		double[] currProfile_plus = comp.getReadProfile_plus(rep.getIndex());
				    		double[] currProfile_minus = comp.getReadProfile_minus(rep.getIndex());
			    			
				    		eventCounter[rep.getIndex()]++;
				    		for(int i=0;i<width;i++){
				    			newModel_plus[rep.getIndex()][i]+=currProfile_plus[offsetLeft+i];
				    			newModel_minus[rep.getIndex()][i]+=currProfile_minus[offsetLeft+i];
				    		}
			    		}
			    		
			    		//Smooth the profile and set a new binding distribution in the BindingModel
			    		if(!logKL[rep.getIndex()].isNaN()){
				    		BindingModel model = rep.getBindingModel();
				    		
				    		for (int i=0;i<width;i++)
				    			newModel[x][i] = newModel_plus[x][i]+newModel_minus[x][i];
				    		
				    		//Smooth the binding model
				    		if(config.getSmoothingBMDuringUpdate()){
				    			if (config.getGaussBMSmooth()){
				    				//Gaussian window smoothing
				    				newModel_plus[x] = StatUtil.gaussianSmoother(newModel_plus[x], config.getBindingModelGaussSmoothParam());
									newModel_minus[x] = StatUtil.gaussianSmoother(newModel_minus[x], config.getBindingModelGaussSmoothParam());
									newModel[x] = StatUtil.gaussianSmoother(newModel[x], config.getBindingModelGaussSmoothParam());
				    			}else if ((int)config.getBindingModelSplineSmoothParam()>0){
						    		// smooth the model profile using spline (not entirely stable)
									newModel_plus[x] = StatUtil.cubicSpline(newModel_plus[x], (int)config.getBindingModelSplineSmoothParam(), (int)config.getBindingModelSplineSmoothParam());
									newModel_minus[x] = StatUtil.cubicSpline(newModel_minus[x], (int)config.getBindingModelSplineSmoothParam(), (int)config.getBindingModelSplineSmoothParam());
									newModel[x] = StatUtil.cubicSpline(newModel[x], (int)config.getBindingModelSplineSmoothParam(), (int)config.getBindingModelSplineSmoothParam());
								}
				    		}
				    		StatUtil.mutate_normalize(newModel_plus[x]);
							StatUtil.mutate_normalize(newModel_minus[x]);
							StatUtil.mutate_normalize(newModel[x]);
				    		//compare models from 2 strands, shift binding position if needed
							//BindingModel.minKL_Shift(newModel_plus[x], newModel_minus[x]);
							
							//use strand-mirrored model for now
							List<Pair<Integer, Double>> dist = new ArrayList<Pair<Integer, Double>>();
							for (int i=0;i<width;i++){
								double sum = newModel[x][i];
								sum = sum>=0?sum:2.0E-300;
								Pair<Integer, Double> p = new Pair<Integer, Double>(i-left, sum);
								dist.add(p);
							}
				
							double[] oldModel = model.getProbabilities();
							model = new BindingModel(dist);
							String outFile = distribFilename+"_"+rep.getName()+".txt";
							model.setFileName(outFile);
							model.printToFile(outFile);
							rep.setBindingModel(model);
							
							logKL[rep.getIndex()] = 0.0;
							if (oldModel.length==width){
								logKL[rep.getIndex()] = StatUtil.log_KL_Divergence(oldModel, model.getProbabilities());
								System.err.println("Refined "+rep.getName()+" read distribution from " + eventCounter[rep.getIndex()] +" binding events.");
							}else{
								System.err.println("Updated "+rep.getName()+" read distribution from " + eventCounter[rep.getIndex()] +" binding events.");
							}
			    		}
	    			}
	    		}
	    	}
	    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions())
	    		cond.updateMaxInfluenceRange();
		}else{
			System.err.println("Read distribution updates turned off.");
		}
		return logKL;
	}
    
    /**
     * Update condition backgrounds for alphas 
     */
    public void updateAlphas(){
    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
			conditionBackgrounds.put(cond, new BackgroundCollection());
			conditionBackgrounds.get(cond).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), cond.getTotalSignalCount(), config.getGenome().getGenomeLength(), config.getMappableGenomeProp(), cond.getMaxModelRange(), '.', 1, true));
			System.err.println("Alpha "+cond.getName()+"\tRange="+cond.getMaxModelRange()+"\t"+(double)conditionBackgrounds.get(cond).getMaxThreshold('.'));
		}
    }
    
    /**
     * Run motif-finding, given the current BindingComponents
     */
    public void updateMotifs(){
    	if(config.getFindingMotifs()){
    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions())
    			motifFinder.findMotifs(cond, activeComponents, trainingRound);
    		motifFinder.alignMotifs();
    		
    		//Print progress
    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
    			if(cond.getFreqMatrix()!=null)
    				System.err.println(cond.getName()+"\t"+WeightMatrix.getConsensus(cond.getFreqMatrix())+"\toffset:"+cond.getMotifOffset());
    			else
    				System.err.println(cond.getName()+"\tNOMOTIF");
    		}
    	}
    }

    /**
     * Initialize the global noise parameters, using only non-potential region counts
     */
    protected void initializeGlobalNoise(){
    	for(int e=0; e<manager.getNumConditions(); e++){
    		ExperimentCondition cond = manager.getExperimentSet().getIndexedCondition(e);
    		double potRegLengthTotal = potRegFilter.getPotRegionLengthTotal();
    		double nonPotRegLengthTotal = config.getGenome().getGenomeLength() - potRegLengthTotal;
    		
    		//Combine control channel counts (avoiding duplication)
    		double potRegCountsSigChannel=0, nonPotRegCountsSigChannel=0; 
    		double potRegCountsCtrlChannel=0, nonPotRegCountsCtrlChannel=0; 
    		ArrayList<Sample> ctrls = new ArrayList<Sample>();
    		for(ControlledExperiment rep : cond.getReplicates())
    			if(rep.hasControl() && !ctrls.contains(rep.getControl())){
    				ctrls.add(rep.getControl());
    				potRegCountsCtrlChannel+=potRegFilter.getPotRegCountsCtrlChannel(rep);
    				nonPotRegCountsCtrlChannel+=potRegFilter.getNonPotRegCountsCtrlChannel(rep);
    			}
    		for(ControlledExperiment rep : cond.getReplicates()){
    			potRegCountsSigChannel+=potRegFilter.getPotRegCountsSigChannel(rep);
				nonPotRegCountsSigChannel+=potRegFilter.getNonPotRegCountsSigChannel(rep);
    		}
    		noisePerBase[e] = nonPotRegCountsSigChannel/nonPotRegLengthTotal;  //Signal channel noise per base
    		System.err.println("Global noise per base initialization for "+cond.getName()+" = "+noisePerBase[e]);
    		//relativeCtrlNoise just tells us if there is a systemic over/under representation of reads in potential regions (in the control)
    		//NOTE: not used for anything right now. 
    		relativeCtrlNoise[e] = (potRegCountsCtrlChannel==0 && nonPotRegCountsCtrlChannel==0) ? 
    				1 : (potRegCountsCtrlChannel/potRegLengthTotal)/(nonPotRegCountsCtrlChannel/nonPotRegLengthTotal);
    	}
    }
    
    /**
     * Update the global noise parameters, using both non-potential region counts and assigned noise responsibilities
     */
    public void updateGlobalNoise(){
    	for(int e=0; e<manager.getNumConditions(); e++){
    		ExperimentCondition cond = manager.getExperimentSet().getIndexedCondition(e);
    		
    		//Don't need to examine noise reads in the update
    		double noiseReads=0; 
    		for(ControlledExperiment rep : cond.getReplicates())
    			noiseReads+=potRegFilter.getNonPotRegCountsSigChannel(rep);

    		for(Region r : noiseResp.keySet())
    			noiseReads+=noiseResp.get(r)[e];
    		noisePerBase[e] = noiseReads/config.getGenome().getGenomeLength();  //Signal channel noise per base
    	}
    }
	
    /**
     * Print all components active at the current time to a file.
     * TESTING ONLY 
     */
    public void printActiveComponentsToFile(){
    	try {
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_t"+trainingRound+".components";
			FileWriter fout = new FileWriter(filename);
			for(Region rr : activeComponents.keySet()){
	    		List<List<BindingComponent>> comps = activeComponents.get(rr);
	    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
	    			for(BindingComponent comp : comps.get(cond.getIndex())){
	    				fout.write(rr.getLocationString()+"\t"+cond.getName()+"\t"+comp.toString()+"\n");			
	    			}
	    		}
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * Print all components active at the current time.
     * TESTING ONLY 
     */
    public void printActiveComponents(){
    	for(Region rr : activeComponents.keySet()){
    		List<List<BindingComponent>> comps = activeComponents.get(rr);
    		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
    			for(BindingComponent comp : comps.get(cond.getIndex())){
    				System.err.println(rr.getLocationString()+"\t"+cond.getName()+"\t"+comp.toString());			
    			}
    		}
    	}
    }
    
	/**
	 * BindingMixtureThread: run binding mixtures over a subset of regions
	 * uniformComponents: if true, components are initialized at uniform prob and spacing, 
	 * 					otherwise, components are initialized from activeComponents
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	class BindingMixtureThread implements Runnable{
		private Collection<Region> regions;
		private int numBindingComponents=1;	//Assuming that the total number of components (active+inactive) is the same in every condition makes coding easier in the BindingEM class.  
		private boolean runEM = true;
		private boolean uniformBindingComponents=false;
		
		public BindingMixtureThread(Collection<Region> regs, boolean EM, boolean uniformBindingComponents){
			regions = regs;	
			this.uniformBindingComponents = uniformBindingComponents;
			runEM=EM;
		}
		
		/**
		 * Run the binding mixture over each test region
		 */
		public void run() {
			//For each region in the test set
        	for (Region rr : regions) {
        		try{
        			//Initialize array of binding component lists, indexed by condition
        			List<List<BindingComponent>> currComps = new ArrayList<List<BindingComponent>>();
        			for(int e=0; e<manager.getNumConditions(); e++)
        				currComps.add(new ArrayList<BindingComponent>());
        			
        			
                    ArrayList<Region> windows = new ArrayList<Region>();
                    if (rr.getWidth()<=config.getBMAnalysisWindowMax())
                        windows.add(rr);
                    else{
                    	//Instead of splitting long windows here, we're splitting in potentialRegionScanner instead.
                    	windows.add(rr);
                    }
                    
                    if(runEM){
		        		//Run EM
	                    Double[] noiseRSums = new Double[manager.getNumConditions()];
	                    for(int e=0; e<manager.getNumConditions(); e++){ noiseRSums[e]=0.0;}
		        		for (Region w : windows){
		        			Pair<List<NoiseComponent>, List<List<BindingComponent>>> wComps = analyzeWindowEM(w);
		        			for(int e=0; e<manager.getNumConditions(); e++){
		        				noiseRSums[e] += wComps.car().get(e).getSumResponsibility();
		        				currComps.get(e).addAll(wComps.cdr().get(e));
		        			}
	                    }
		        		
		        		//Only non-zero components are returned by analyzeWindow, so add them to the recorded active components
		        		synchronized(activeComponents){	activeComponents.put(rr, currComps);}
		        		
		        		//Add the sum of noise responsibilities to this region
		        		synchronized(noiseResp){ noiseResp.put(rr, noiseRSums);}
                    }else{
                    	//Run ML assignment
                    	List<BindingEvent> windowBindingEvents = new ArrayList<BindingEvent>();
                    	for (Region w : windows){
                    		windowBindingEvents.addAll(  analyzeWindowML(w) );
                    	}
                    	synchronized(bindingEvents){bindingEvents.addAll(windowBindingEvents);}
                    }
        		} catch(Exception e){
                    System.err.println("ERROR: Exception when analyzing region "+rr.toString());
                    e.printStackTrace(System.err);
                    System.exit(-1);
                }
            }
		}
		
		/**
		 * Train a BindingComponent EM over a given window
		 * 
		 * In the old implementation, there was a choice here to run EM or scan for the ML soln
         * (if we could be sure there was only one binding event). However, we can no longer scan because
         * the solution can depend on multiple conditions. So we always do EM.
         * 
         * In addition, since we are now using components that can change their position along the genome, 
         * there is no need for the component resolution updates. 
         * 
         * Components are now semi-independently estimated for each condition, so they have become lists-of-lists.
         * 
         * We also now initialize noise components, which are position-less and have a fixed (estimated) emission probability.
         *  
		 * @param w
		 * @return Pair of component lists (noise components and binding components) indexed by condition
		 */
		private Pair<List<NoiseComponent>, List<List<BindingComponent>>> analyzeWindowEM(Region w){
			BindingEM EM = new BindingEM(config, manager, conditionBackgrounds, potRegFilter.getPotentialRegions().size());
			List<List<BindingComponent>> bindingComponents=null;
			List<NoiseComponent> noiseComponents=null;
			List<List<BindingComponent>> nonZeroComponents = new ArrayList<List<BindingComponent>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				nonZeroComponents.add(new ArrayList<BindingComponent>());
			
			//Check if this is a region to plot (and if so, calculate any reqd offsets)
			Region plotSubReg=null;
			for(Region p : regionsToPlot)
				if(w.overlaps(p))
					plotSubReg = p;
			
			//Load signal data
			List<List<StrandedBaseCount>> signals = loadSignalData(w);
            if (signals==null)
                return new Pair<List<NoiseComponent>, List<List<BindingComponent>>>(noiseComponents, nonZeroComponents);
            //Load control data
            List<List<StrandedBaseCount>> controls = loadControlData(w);
            
            //Initialize noise components
            noiseComponents = initializeNoiseComponents(w, signals, controls);

            //Initialize binding components
            if(uniformBindingComponents)
            	bindingComponents = initializeBindingComponentsUniformly(w, noiseComponents);
            else
            	bindingComponents = initializeBindingComponentsFromAllConditionActive(w, noiseComponents, true);
            
            //Motif prior
            String seq = config.getFindingMotifs() ? motifFinder.getSeq(w):null;
            double[][] motifPrior = config.getFindingMotifs() ? motifFinder.scanRegionWithMotifs(w, seq) : null;
            
            //EM learning: resulting binding components list will only contain non-zero components
            nonZeroComponents = EM.train(signals, w, noiseComponents, bindingComponents, numBindingComponents, motifPrior, trainingRound, plotSubReg);
           
            return new Pair<List<NoiseComponent>, List<List<BindingComponent>>>(noiseComponents, nonZeroComponents);
        }//end of analyzeWindowEM method
		
		/**
		 * Assign BindingComponents over a given window with ML solution
		 *  
		 * @param w
		 * @return Pair of component lists (noise components and binding components) indexed by condition
		 */
		private List<BindingEvent> analyzeWindowML(Region w){
			BindingMLAssignment ML = new BindingMLAssignment(config, manager, conditionBackgrounds, potRegFilter.getPotentialRegions().size());
			List<BindingComponent> bindingComponents=null;
			List<NoiseComponent> noiseComponents=null;
			List<BindingEvent> currEvents = new ArrayList<BindingEvent>(); 
			
			//Load signal data
			List<List<StrandedBaseCount>> signals = loadSignalData(w);
            if (signals==null)
                return currEvents;
            //Load control data
            List<List<StrandedBaseCount>> controls = loadControlData(w);
            
            //Initialize noise components
            noiseComponents = initializeNoiseComponents(w, signals, controls);

            
            //Configurations seen in another condition
            ArrayList<ComponentConfiguration> seenConfigs = new ArrayList<ComponentConfiguration>();
    		//Assign reads to components
            for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
            	//Initialize binding components: shared configuration or condition-specific
            	if(config.getMLSharedComponentConfiguration()){
            		bindingComponents = initializeBindingComponentsFromAllConditionActive(w, noiseComponents, false).get(cond.getIndex());
            	}else{
            		bindingComponents = initializeBindingComponentsFromOneConditionActive(w, noiseComponents.get(cond.getIndex()), cond.getIndex());
            	}
            	int numComp = bindingComponents.size();
            	
            	//Construct configuration
    			ComponentConfiguration currCC = new ComponentConfiguration(bindingComponents, cond.getIndex());
    		
    			//Have we already seen this configuration?
    			boolean ccFound=false;
    			for(ComponentConfiguration cc : seenConfigs){
    				if(currCC.isSameAs(cc)){//If so, add another valid condition to the observed binding events
    					int parent = cc.getParentCondition();
    					if(!config.getMLSharedComponentConfiguration()){
    						for(BindingEvent be : currEvents)
    							if(be.isFoundInCondition(parent))
    								be.setIsFoundInCondition(cond.getIndex(),true);
    					}
    					ccFound=true;
    					break;
    				}
    			}
    			
    			if(!ccFound){
    				//Add configuration to seen
    				seenConfigs.add(currCC);
    				
    				//ML assignment
    				List<BindingEvent> condEvents = ML.assign(signals, controls, w, noiseComponents, bindingComponents, numComp);
    				for(BindingEvent be : condEvents)
    					if(config.getMLSharedComponentConfiguration())
    						setFoundInConditions(be, w);
    					else
    						be.setIsFoundInCondition(cond.getIndex(),true);
    				currEvents.addAll(condEvents);
    			}
            }
            
            //If we haven't used shared component ML, we need to edit and consolidate binding events
            // 1) Consolidate, because otherwise you can have duplicate binding events
            // 2) Edit - set counts to zero at conditions where the event is not active
            if(!config.getMLSharedComponentConfiguration()){
            	currEvents = consolidateBindingEvents(currEvents);
            }
            
            //Add in sequences and final motif scores here
            if(config.isAddingSequences()){
	            String seq = config.getFindingMotifs() ? motifFinder.getSeq(w):null;
	            Pair<Double[][], String[][]> motifScores = config.getFindingMotifs() ? motifFinder.scanRegionWithMotifsGetSeqs(w, seq) : null;
	            if(seq!=null && motifScores!=null){
		            for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
		            	for(BindingEvent b: currEvents){
		            		b.setMotifScore(cond, motifScores.car()[cond.getIndex()][b.getPoint().getLocation()-w.getStart()]);
		            		b.setSequence(cond, motifScores.cdr()[cond.getIndex()][b.getPoint().getLocation()-w.getStart()]);
		            	}
		            }
	            }
            }
            
            return currEvents;
        }//end of analyzeWindowML method

		/**
		 * Set the conditions in which a given binding event is still active after EM training
		 * @param b
		 * @param currReg
		 */
		private void setFoundInConditions(BindingEvent b, Region currReg){
			for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingComponent comp : activeComponents.get(currReg).get(e)){
        			if(comp.getPosition() == b.getPoint().getLocation()){
        				b.setIsFoundInCondition(e, true);
        			}
        		}
		}
		
		/**
		 * Load all signal read hits in a region by condition. 
		 * 
		 * @param w
		 * @return List of List of StrandedBaseCounts, indexed by replicate index
		 */
		private List<List<StrandedBaseCount>> loadSignalData(Region w){
			List<List<StrandedBaseCount>> data = new ArrayList<List<StrandedBaseCount>>();
			for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
				for(ControlledExperiment rep : cond.getReplicates())
					data.add(new ArrayList<StrandedBaseCount>());
			}
			for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
				for(ControlledExperiment rep : cond.getReplicates()){
					data.get(rep.getIndex()).addAll(rep.getSignal().getUnstrandedBases(w));
				}
			}
			return data;
		}
		
		/**
		 * Load all control read hits in a region by condition. 
		 * 
		 * @param w
		 * @return List of List of StrandedBaseCounts, indexed by replicate index
		 */
		private List<List<StrandedBaseCount>> loadControlData(Region w){
			List<List<StrandedBaseCount>> data = new ArrayList<List<StrandedBaseCount>>();
			for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
				for(ControlledExperiment rep : cond.getReplicates())
					data.add(new ArrayList<StrandedBaseCount>());
			}
			for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
				for(ControlledExperiment rep : cond.getReplicates()){
					if(rep.hasControl())
						data.get(rep.getIndex()).addAll(rep.getControl().getUnstrandedBases(w));
				}
			}
			return data;
		}
		
		/**
         * Initializes the components uniformly: i.e. space them evenly along the region.
         *
         * @param currReg
         */
        private List<List<BindingComponent>> initializeBindingComponentsUniformly(Region currReg, List<NoiseComponent> noise){
        	List<List<BindingComponent>> components = new ArrayList<List<BindingComponent>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				components.add(new ArrayList<BindingComponent>());
			
            //Place components along region
			int componentSpacing =  config.INIT_COMPONENT_SPACING;
			if(componentSpacing >= currReg.getWidth())
				System.err.println("Error:  region width less than component spacing in "+currReg.getLocationString());

			numBindingComponents = currReg.getWidth()/componentSpacing;

            //Set up the components
			for(int e=0; e<manager.getNumConditions(); e++){
	            double numC=0;
	            for(int i=0; i<numBindingComponents; i++){
	                Point pos = new Point(config.getGenome(), currReg.getChrom(), currReg.getStart()+(i*componentSpacing));
	                BindingComponent currComp = new BindingComponent(pos, manager.getExperimentSet().getReplicates().size());
	                currComp.setIndex(i);
	                numC++;
	                components.get(e).add(currComp);
	            }
	            //Initialize normalized mixing probabilities (subtracting the noise emission probability)
	            double emission = (1-noise.get(e).getPi())/numC;
	            for(BindingComponent b : components.get(e)){
	                b.uniformInit(emission);	                
	            }
			}
            return components; 
        }//end of initializeComponents method

        /**
         * Initializes the components from all active components in all conditions: 
         * 		Uses all active component locations from the last round of training in each condition,
         * 		i.e. not just the conditions in which those active components were active. Also adds in 
         * 		extra components flanking the active locations in case the binding distribution update
         * 		has made more joint events separable. If no components exist, a rescue component is added. 
         *
         * @param currReg
         */
        private List<List<BindingComponent>> initializeBindingComponentsFromAllConditionActive(Region currReg, List<NoiseComponent> noise, boolean addFlanking){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingComponent comp : activeComponents.get(currReg).get(e)){
        			if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd())
        				componentPositions.add(comp.getPosition());
        			if(addFlanking){
        				if(!componentPositions.contains(comp.getPosition()-config.getAddFlankingComponentSpacing())
        						 && comp.getPosition()-config.getAddFlankingComponentSpacing()>=currReg.getStart())
        					componentPositions.add(comp.getPosition()-config.getAddFlankingComponentSpacing());
        				if(!componentPositions.contains(comp.getPosition()+config.getAddFlankingComponentSpacing())
        						&& comp.getPosition()+config.getAddFlankingComponentSpacing()<currReg.getEnd())
        					componentPositions.add(comp.getPosition()+config.getAddFlankingComponentSpacing());
        			}
        		}

        	numBindingComponents = componentPositions.size();
        	
        	//If no components exist in region, add one to the center to allow rescues
        	if(numBindingComponents==0 && addFlanking){
        		componentPositions.add(currReg.getMidpoint().getLocation());
        		numBindingComponents++;
        	}

        	//Make new components with these locations
        	List<List<BindingComponent>> components = new ArrayList<List<BindingComponent>>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		components.add(new ArrayList<BindingComponent>());
        	
        	//Set up the components
        	for(int e=0; e<manager.getNumConditions(); e++){
        		double numC=(double)numBindingComponents; int index=0;
        		double emission = (1-noise.get(e).getPi())/numC;
	    		for(Integer i : componentPositions){
	    			Point pos = new Point(config.getGenome(), currReg.getChrom(), i);
	    			BindingComponent currComp = new BindingComponent(pos, manager.getExperimentSet().getReplicates().size());
	    			currComp.setIndex(index);
	    			index++;
	    			
    				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
        			currComp.uniformInit(emission);
    				components.get(e).add(currComp);
    			}
    		}
        	return components; 
        }//end of initializeComponents method

        /**
         * Initializes components from active components in a single condition: 
         * 		Uses active component locations from the last round of training in one condition,
         * 		No flanking components or resuce components added here, since resulting components will only be used
         * 		in ML assignment.  
         *
         * @param currReg
         */
        private List<BindingComponent> initializeBindingComponentsFromOneConditionActive(Region currReg, NoiseComponent noise, int conditionIndex){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(BindingComponent comp : activeComponents.get(currReg).get(conditionIndex)){
        		if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd())
        			componentPositions.add(comp.getPosition());
        	}

        	numBindingComponents = componentPositions.size();

        	//Make new components with these locations
        	List<BindingComponent> components = new ArrayList<BindingComponent>();
        	
        	//Set up the components
        	double numC=(double)numBindingComponents; int index=0;
    		double emission = (1-noise.getPi())/numC;
    		for(Integer i : componentPositions){
    			Point pos = new Point(config.getGenome(), currReg.getChrom(), i);
    			BindingComponent currComp = new BindingComponent(pos, manager.getExperimentSet().getReplicates().size());
    			currComp.setIndex(index);
    			index++;
    			//Initialize normalized mixing probabilities (subtracting the noise emission probability)
    			currComp.uniformInit(emission);
				components.add(currComp);
			}
    		
        	return components; 
        }//end of initializeComponents method
        
        /**
         * Initializes the noise components.
         * Noise distributions are set per replicate, so control channel reads are not combined into a single list. 
         *
         * @param currReg
         */
        private List<NoiseComponent> initializeNoiseComponents(Region currReg, List<List<StrandedBaseCount>> sigHits, List<List<StrandedBaseCount>> ctrlHits){
        	List<NoiseComponent> noise = new ArrayList<NoiseComponent>();
        	int numReps = manager.getExperimentSet().getReplicates().size();
        	double [] localSigRepCounts=new double [numReps];
        	double [] localCtrlRepCounts=new double [numReps];
        	
        	//Calculate expected noise distributions
        	double [][] distribs=new double[numReps][];
        	for(ExperimentCondition cond : manager.getExperimentSet().getConditions())
        		for(ControlledExperiment rep : cond.getReplicates()){
	    			if(rep.hasControl() && ctrlHits.get(rep.getIndex()).size()>0){
	            		distribs[rep.getIndex()] = smoothNoiseDistribs(currReg, ctrlHits.get(rep.getIndex()));
	            		localCtrlRepCounts[rep.getIndex()]=0;
	            		for(StrandedBaseCount b : ctrlHits.get(rep.getIndex()))
	            			localCtrlRepCounts[rep.getIndex()]+=b.getCount();
	    			}else
	    				distribs[rep.getIndex()] = null;
	    		}
        	
        	//Initialize the noise component
        	for(int e=0; e<manager.getNumConditions(); e++){
        		ExperimentCondition cond = manager.getExperimentSet().getIndexedCondition(e);
        		double emission = 0;
        		
        		//Sum signal reads & set local region experiment counts
        		double sigCounts=0;
        		for(ControlledExperiment rep : cond.getReplicates()){
        			localSigRepCounts[rep.getIndex()]=0;
        			for(StrandedBaseCount b : sigHits.get(rep.getIndex())){
        				sigCounts+=b.getCount(); localSigRepCounts[rep.getIndex()]+=b.getCount();
        			}
        		}
        		
        		//Calculate a local noise factor to check for expected over-representation of noise reads, as specified in the control.
        		//We assume that only noise generates the control channel. Therefore, the first two terms in the localNoiseFactor calculation
        		//test for over-representation in the observed control read counts in the local window. 
        		//The last term in the calculation is a local replicate weight, to account for the fact that some replicate signal channels are 
        		//producing more of the signal reads (and of course, the noise we are actually accounting for is in the signal channel). 
        		double localNoiseFactor = 0;
        		for(ControlledExperiment rep : cond.getReplicates()){
        			if(rep.hasControl() && ctrlHits.get(rep.getIndex()).size()>0){
        				localNoiseFactor+=(localCtrlRepCounts[rep.getIndex()]/(double)currReg.getWidth()) /
        								  (rep.getControl().getHitCount()/(double)currReg.getWidth())     *
        								  (localSigRepCounts[rep.getIndex()]/sigCounts); //over-rep x weight
        			}else{
        				localNoiseFactor+=(localSigRepCounts[rep.getIndex()]/sigCounts); //1 x weight
        			}
        		}
        		
        		//Calculate expected noise emission = number of expected noise reads over the total signal reads in this region
        		// If local noise factor is above 1, account for it. Otherwise, meh.
        		emission = (noisePerBase[e] * (double)currReg.getWidth()) / sigCounts;
        		if(localNoiseFactor>1)
        			emission*=localNoiseFactor;
        		if(emission>config.NOISE_EMISSION_MAX)
        			emission = config.NOISE_EMISSION_MAX;
        		if(emission<config.NOISE_EMISSION_MIN)
        			emission = config.NOISE_EMISSION_MIN;
        		
        		//Add the noise component
        		NoiseComponent n = new NoiseComponent(emission, distribs, currReg, numReps);
        		noise.add(n);
        	}
        	return noise;
        }//end of initializeNoiseComponents method
      
        /**
         * Smooth the distribution of the control reads over the window to make a probability distribution of noise over the region
         * @param currReg
         * @param ctrlHits
         * @return double array the same width as the region, containing probabilities normalized to sum to 1
         */
        private double[] smoothNoiseDistribs(Region currReg, List<StrandedBaseCount> ctrlHits){
        	double [] distrib = new double[currReg.getWidth()];
        	double [] counts = new double[currReg.getWidth()];
        	//Pseudocounts for distrib
        	for(int d=0; d<currReg.getWidth(); d++)
        		counts[d]=1;
        	//Add in count weights
        	for(StrandedBaseCount hit : ctrlHits){
        		int index = hit.getCoordinate()-currReg.getStart();
        		if(index>=0 && index<currReg.getWidth())
        			counts[index]+=hit.getCount();
        	}
        	
        	//Smooth
        	for(int i=0; i<currReg.getWidth(); i++){
        		double sum=0, num=0;
        		for(int s=i-(config.NOISE_DISTRIB_SMOOTHING_WIN/2); s<i+(config.NOISE_DISTRIB_SMOOTHING_WIN/2); s++){
        			if(s>=0 && s<currReg.getWidth()){
        				num++; sum+=counts[s];
        			}
        		}
        		distrib[i] = sum/num;
        	}
        	
        	//Normalize
        	double total = 0;
        	for(int d=0; d<currReg.getWidth(); d++)
        		total+=distrib[d];
        	for(int d=0; d<currReg.getWidth(); d++)
        		distrib[d] = distrib[d]/total;
        	
        	return distrib;
        }
	}
	
	/**
	 * Consolidate and edit binding events.
	 * Follows a strict definition of binding event quantification - 
	 * if the event is not present in the condition, it gets a zero count assigned.
	 * Also merges positional duplicate events. 
	 * @param ev
	 * @return
	 */
	private List<BindingEvent> consolidateBindingEvents(List<BindingEvent> ev){
		List<BindingEvent> newEvents = new ArrayList<BindingEvent>();
		HashMap<Point, Integer> eventMap = new HashMap<Point, Integer>();
		int count=0;
		for(BindingEvent be : ev){
			if(!eventMap.containsKey(be.getPoint())){
				eventMap.put(be.getPoint(), count);
				count++;
				
				newEvents.add(be);
				//First time an event is added, clear out inactive events
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					if(!be.isFoundInCondition(cond)){
						be.setCondSigHits(cond, 0.0);
		        		for(ControlledExperiment rep : cond.getReplicates()){
		        			be.setRepSigHits(rep, 0.0);
		        		}if(config.CALC_COMP_LL)
			            	be.setLLd(cond, 0);
		        }	}
			}else{
				int index = eventMap.get(be.getPoint());
				//For events that are already in the list, just update the active condition read counts
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					if(be.isFoundInCondition(cond)){
						newEvents.get(index).setCondSigHits(cond, be.getCondSigHits(cond));
						newEvents.get(index).setCondCtrlHits(cond, be.getCondCtrlHits(cond));
		        		for(ControlledExperiment rep : cond.getReplicates()){
		        			newEvents.get(index).setRepSigHits(rep, be.getRepSigHits(rep));
							newEvents.get(index).setRepCtrlHits(rep, be.getRepCtrlHits(rep));
		        		}if(config.CALC_COMP_LL)
		        			newEvents.get(index).setLLd(cond, be.getLLd(cond));
		        }	}
			}
		}
		return newEvents;
	}
	
	/**
	 * ComponentConfiguration: represents a configuration of binding components as an array of positions.
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	protected class ComponentConfiguration{
		int [] positions=null;
		int parentCondIndex;
		//Constructor
		public ComponentConfiguration(List<BindingComponent> comps, int parentCondition){
			Collections.sort(comps);
			positions = new int[comps.size()];
			for(int p=0; p<comps.size(); p++)
				positions[p]=comps.get(p).getPosition();
			parentCondIndex = parentCondition;
		}
		//Return the index of the condition that initialized this configuration
		public int getParentCondition(){return parentCondIndex;}
		//Compare two configurations
		public boolean isSameAs(ComponentConfiguration cc){
			if(positions.length != cc.positions.length)
				return false;
			boolean isEqual=true;
			for(int p=0; p<positions.length; p++)
				isEqual = isEqual && positions[p]==cc.positions[p];
			return isEqual;
		}
	}
}
