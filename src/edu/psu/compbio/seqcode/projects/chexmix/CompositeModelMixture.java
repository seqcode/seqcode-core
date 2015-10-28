package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.fitting.GaussianFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;

/**
 * CompositeModelMixture: defines a protein-DNA interaction mixture model in a composite tag distribution.
 * 
 *  This mixture assumes that the same protein-DNA interaction model is valid across all examined conditions. 
 * 
 * TODO: support multi-rep TagDistributions
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CompositeModelMixture {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ChExMixConfig config;
	protected ExperimentManager manager;
	protected CompositeTagDistribution compositeDistrib; //The composite tag distribution under investigation
	protected CompositeTagDistribution controlCompositeDistrib; //Optional composite tag distribution from matching control
	protected ProteinDNAInteractionModel model; //The model to train
	protected List<CompositeModelSiteAssignment> siteAssignments;
	protected CompositeModelEM EMtrainer; //Training method
	
	protected TagDistribution initBackDistrib, initCSDistrib, initXODistrib;
	protected double initNoisePi;
	protected List<Region> regionsToPlot;
	protected int trainingRound=0;
	
	public CompositeModelMixture(CompositeTagDistribution tagDist, CompositeTagDistribution ctrlDist, GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon, ExperimentManager eMan){
		gconfig = gcon;
		econfig = econ;
		config = ccon;
		manager = eMan;
		compositeDistrib = tagDist;
		controlCompositeDistrib = ctrlDist;
		regionsToPlot = config.getRegionsToPlot();

		EMtrainer = new CompositeModelEM(compositeDistrib, config, manager);

		initializeTagDistributions();
		
		model = new ProteinDNAInteractionModel(config, compositeDistrib.getWinSize(), initXODistrib, initCSDistrib, initBackDistrib, initNoisePi);
		siteAssignments = new ArrayList<CompositeModelSiteAssignment>();
	}
	
	//Accessors
	public ProteinDNAInteractionModel getModel(){return model;}
	
	/**
	 * Run the mixture model  
	 * 
	 * EM: if true, run EM on composite, otherwise run ML assignment on component points
	 */
	public void execute(boolean EM){
		trainingRound++;
		
		try{
			if(EM){
				//Run EM outer loops
				boolean converged=false;
				double[] kl = new double[config.getMaxModelUpdateRounds()];
				for(trainingRound=0; trainingRound<config.getMaxModelUpdateRounds() && !converged; trainingRound++){
					//Run EM inner loops
					model = EMtrainer.train(model, trainingRound);
					kl[trainingRound]=updateTagDistributions();
				
					System.out.println("TrainingRound: "+trainingRound+":\n"+model.toString()+"\n\tModelKL="+kl[trainingRound]);
					
					//Check for convergence
					if(kl[trainingRound]<config.getModelConvergenceKL())
						converged=true;
				}
			}else{
				//TODO: multithreading
				
				//List all site indexes (looks dumb now, but this will enable multithreading later)
				List<Integer> allSites = new ArrayList<Integer>();
				for(int s=0; s<compositeDistrib.getPoints().size(); s++)
					allSites.add(s);
				
				//Run ML
				CompositeModelMLAssign MLassigner = new CompositeModelMLAssign(compositeDistrib, allSites, config, manager);
				List<CompositeModelSiteAssignment> currAssignments = MLassigner.assign(model);
				
				synchronized(siteAssignments){
					siteAssignments.addAll(currAssignments);
					Collections.sort(siteAssignments);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
     
	/**
	 * Initialize the CS, XO, and background tag distributions
	 */
	protected void initializeTagDistributions(){
		//CS (empirical)
		/*if(config.getDefaultCSModel()!=null)
			initCSDistrib = config.getDefaultCSModel();
		else
			initCSDistrib = new TagDistribution(TagDistribution.defaultChipSeqEmpiricalDistribution, null);
		*/
		initCSDistrib = new TagDistribution(compositeDistrib.getWinSize());
		initCSDistrib.loadGaussianDistrib(150, 100);
		
		//XO
		initXODistrib = new TagDistribution(200);
		initXODistrib.loadGaussianDistrib(config.getXLDistribOffset(), config.getXLDistribSigma());
		
		//Background
		if(controlCompositeDistrib==null){
			initBackDistrib = new TagDistribution(compositeDistrib.getWinSize());
			initBackDistrib.loadFlatDistrib();
			//estimate noise rate from the tails of the composite distribution
			double sum=0, count=0;
			for(ExperimentCondition cond : manager.getConditions()){
				double[] wprobs=compositeDistrib.getCompositeWatson(cond), cprobs=compositeDistrib.getCompositeCrick(cond);
				for(int x=0; x<10; x++){
					sum+=wprobs[wprobs.length-x-1];
					sum+=cprobs[x];
					count+=2;
				}
			}
			initNoisePi=Math.min(config.NOISE_EMISSION_MAX, Math.max(config.NOISE_EMISSION_MIN, (sum/count)*(compositeDistrib.getWinSize())));
		}else{
			//Make a paired list from the summed counts in the control tag distribution
			double[] wcounts = new double[controlCompositeDistrib.getWinSize()];
			double[] ccounts = new double[controlCompositeDistrib.getWinSize()];
			double sum=0;
			for(ExperimentCondition cond : manager.getConditions()){
				double[] wtmp = controlCompositeDistrib.getCompositeWatson(cond);
				double[] ctmp = controlCompositeDistrib.getCompositeWatson(cond);
				for(int x=0; x<controlCompositeDistrib.getWinSize(); x++){
					wcounts[x]+=wtmp[x];
					ccounts[x]+=ctmp[x];
					sum+=wtmp[x]+ctmp[x];
				}
			}
			initBackDistrib = new TagDistribution(compositeDistrib.getWinSize());
			try {
				initBackDistrib.loadData(wcounts, ccounts);
			} catch (Exception e) {
				e.printStackTrace();
			}
			initNoisePi=sum/(controlCompositeDistrib.getWinSize()*2);
		}
	}
	
	/**
     * Update tag distributions given the current component responsibilities. 
     * 
     * @return log KL values
     */
    public Double updateTagDistributions(){
    	double logKL = 0;
    	TagDistribution CSdistrib = model.getCSTagDistribution();
    	TagDistribution XLdistrib = model.getXLTagDistribution();
    	
    	//CS component (empirical or gaussian fit) 
    	double[] oldCSModelW = CSdistrib.getWatsonProbabilities();
    	double[] oldCSModelC = CSdistrib.getCrickProbabilities();
		double[] csW = model.getCSComponent().getTagProfile(true);
    	double[] csC = model.getCSComponent().getTagProfile(false);
    	//smooth with arbitrary gaussian
    	csW = StatUtil.gaussianSmoother(csW, 10);
    	csC = StatUtil.gaussianSmoother(csC, 10);
    	/*
    	double[] newCSModelW=new double[CSdistrib.getWinSize()];
		double[] newCSModelC=new double[CSdistrib.getWinSize()];
		for(int i=0; i<CSdistrib.getWinSize(); i++){
			newCSModelW[i]=0; newCSModelC[i]=0;
		}int i=0;
    	for(int x=offsetLeft; x<=offsetRight; x++){
    		newCSModelW[i]=csW[x];
    		newCSModelC[i]=csC[x];
    		i++;
    	}
    	StatUtil.mutate_normalize(newCSModelW);
    	StatUtil.mutate_normalize(newCSModelC);
    	try {
			CSdistrib.loadData(newCSModelW, newCSModelC);
		} catch (Exception e) {
			e.printStackTrace();
		}*/
    	GaussianFitter fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
    	for(int i=0; i<CSdistrib.getWinSize(); i++){
    		fitter.addObservedPoint((double)(i+CSdistrib.getLeft()), csW[i]+csC[CSdistrib.getWinSize()-i-1]);
    	}
    	double[] parameters = fitter.fit();;
    	double newOffset = -1*parameters[1];
    	double newSigma = parameters[2];
    	System.out.println("CSGaussianFit:\t"+newOffset+"\t"+newSigma);
    	CSdistrib.loadGaussianDistrib(newOffset, newSigma);
    	double[] newCSModelW = CSdistrib.getWatsonProbabilities();
    	double[] newCSModelC = CSdistrib.getCrickProbabilities();
    	
    	
    	//XL components (fit gaussian)
    	double[] oldXLModelW = XLdistrib.getWatsonProbabilities();
    	double[] oldXLModelC = XLdistrib.getCrickProbabilities();
		double[] xlW = new double[XLdistrib.getWinSize()];
    	double[] xlC = new double[XLdistrib.getWinSize()];
    	for(int i=0; i<XLdistrib.getWinSize(); i++){ xlW[i]=0; xlC[i]=0;}
    	for(CompositeModelComponent xlComp : model.getXLComponents()){
    		double[] currW = xlComp.getTagProfile(true);
    		double[] currC = xlComp.getTagProfile(false);
    		for(int i=0; i<XLdistrib.getWinSize(); i++){ 
    			xlW[i]+=currW[i]; 
    			xlC[i]+=currC[i];
    		}
    	}
    	//smooth with arbitrary gaussian
    	xlW = StatUtil.gaussianSmoother(xlW, 1);
    	xlC = StatUtil.gaussianSmoother(xlC, 1);
    	fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
    	for(int i=0; i<XLdistrib.getWinSize(); i++){
    		fitter.addObservedPoint((double)(i+XLdistrib.getLeft()), xlW[i]+xlC[XLdistrib.getWinSize()-i-1]);
    	}
    	parameters = fitter.fit();;
    	newOffset = -1*parameters[1];
    	newSigma = parameters[2];
    	System.out.println("XLGaussianFit:\t"+newOffset+"\t"+newSigma);
    	XLdistrib.loadGaussianDistrib(newOffset, newSigma);
    	double[] newXLModelW = XLdistrib.getWatsonProbabilities();
    	double[] newXLModelC = XLdistrib.getCrickProbabilities();
    	
    	
    	//Calc KL
    	logKL += StatUtil.log_KL_Divergence(oldCSModelW, newCSModelW);
    	logKL += StatUtil.log_KL_Divergence(oldCSModelC, newCSModelC);
    	logKL += StatUtil.log_KL_Divergence(oldXLModelW, newXLModelW);
    	logKL += StatUtil.log_KL_Divergence(oldXLModelC, newXLModelC);
    	
		return logKL;
	}
    
}
