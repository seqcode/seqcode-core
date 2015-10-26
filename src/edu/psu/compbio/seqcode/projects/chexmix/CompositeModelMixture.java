package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.List;

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
				model = EMtrainer.train(model, trainingRound);
			}else{
				//runML();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
     
	/**
	 * Initialize the CS, XO, and background tag distributions
	 */
	protected void initializeTagDistributions(){
		//CS
		if(config.getDefaultCSModel()!=null)
			initCSDistrib = config.getDefaultCSModel();
		else
			initCSDistrib = new TagDistribution(TagDistribution.defaultChipSeqEmpiricalDistribution, null);
		
		//XO
		initXODistrib = new TagDistribution(200);
		initXODistrib.loadGaussianDistrib(config.getXLDistribOffset(), config.getXLDistribSigma());
		
		//Background
		if(controlCompositeDistrib==null){
			initBackDistrib = new TagDistribution(compositeDistrib.getWinSize());
			initBackDistrib.loadFlatDistrib();
			//estimate noise rate from the tails of the composite distribution
			double[] wprobs=initBackDistrib.getWatsonProbabilities(), cprobs=initBackDistrib.getCrickProbabilities();
			double sum=0, count=0;
			for(int x=0; x<10; x++){
				sum+=wprobs[wprobs.length-x-1];
				sum+=cprobs[x];
				count+=2;
			}
			initNoisePi=sum/count;
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
    	
    	//CS component (empirical) 
    	double[] oldCSModelW = CSdistrib.getWatsonProbabilities();
    	double[] oldCSModelC = CSdistrib.getCrickProbabilities();
    	int tagProfileCenter = config.MAX_BINDINGMODEL_WIDTH/2;
    	int offsetLeft = CSdistrib.getLeft()+tagProfileCenter;
    	int offsetRight = tagProfileCenter+CSdistrib.getRight();
		double[] newCSModelW=new double[CSdistrib.getWinSize()];
		double[] newCSModelC=new double[CSdistrib.getWinSize()];
		for(int i=0; i<CSdistrib.getWinSize(); i++){
			newCSModelW[i]=0; newCSModelC[i]=0;
		}
		double[] csW = model.getCSComponent().getTagProfile(true);
    	double[] csC = model.getCSComponent().getTagProfile(false);
    	//smooth with arbitrary gaussian
    	csW = StatUtil.gaussianSmoother(csW, 2);
    	csC = StatUtil.gaussianSmoother(csC, 2);
    	int i=0;
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
		}

    	//XL components (fir gaussian)
    	//TODO
    	
    	logKL += StatUtil.log_KL_Divergence(oldCSModelW, newCSModelW);
    	logKL += StatUtil.log_KL_Divergence(oldCSModelC, newCSModelC);
    	
		return logKL;
	}
}
