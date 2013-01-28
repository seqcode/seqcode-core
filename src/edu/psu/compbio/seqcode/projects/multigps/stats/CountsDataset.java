package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Jama.Matrix;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.utilities.ScatterPlotMaker;

/**
 * CountsDataset is a class used to hold a matrix of counts associated with a set of experiments.
 * The class is used as the basis for pairwise statistical testing. As such, the class also 
 * maintains scaling factors (for each pair of experiments), and p-value matrices (for each pair
 * of conditions). 
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CountsDataset {

	protected Matrix counts; //Counts matrix. Indexed by unit and sample
	protected double [] totals;   //Count totals. Indexed by sample
	protected double [] scaling;  //Scaling factors. Indexed by sample
	protected String [] unitNames; //Unit names (gene names or peak coords)
	protected int numSamples=1;    
	protected int numUnits = 1;
	protected int numConds = 1;
	protected int focalCondition=0;
	protected int [] design;		//Design array. Sample index --> Condition index
	protected Matrix DEpval;	//P-values of differential expression. Indexed by unit and condition
	protected Matrix condFold;//Condition Fold (logFold from EdgeR). Indexed by unit and condition
	protected Matrix condMean;//Condition means (logCPM from EdgeR). Indexed by unit and condition
	protected Matrix condRawVar; //Condition variances. Indexed by unit and condition
	protected HashMap<Integer, Pair<String,String>> sampleToExptName;
	protected HashMap<Integer, String> condToName;
	protected HashMap<String, Integer> unitToIndex;
	protected final double LOG_2 = Math.log(2.0);
	
	/**
	 * Constructor: preformed data
	 * @param counts : count matrix
	 * @param units : Names of each unit (i.e. gene name or peak coordinate) 
	 * @param design : design array (which index in counts corresponds to which condition, indexed by integer)
	 * @param exptNameMap: translation between sample index and condition/replicate names
	 * @param condNameMap: translation between condition index and name
	 * @param focalCond: focal condition (sample index of)
	 */
	public CountsDataset(Matrix counts, String [] units, int [] design, HashMap<Integer, Pair<String,String>> exptNameMap, HashMap<Integer, String> condNameMap){this(counts, units, design, exptNameMap, condNameMap, 0);}
	public CountsDataset(Matrix counts, String [] units, int [] design, HashMap<Integer, Pair<String,String>> exptNameMap, HashMap<Integer, String> condNameMap, int focalCond){
		this.counts = counts;
		this.unitNames = units;
		this.design = design;
		this.focalCondition = focalCond;
		this.sampleToExptName = exptNameMap;
		this.numSamples = counts.getColumnDimension();
		this.numUnits = counts.getRowDimension();
		this.condToName = condNameMap;
		this.numConds = condToName.size();
		this.totals = new double[numSamples];
		this.scaling = new double[numSamples];
		this.DEpval = new Matrix(numUnits,numConds);
		this.condFold = new Matrix(numUnits,numConds);
		this.condMean = new Matrix(numUnits,numConds);
		this.condRawVar = new Matrix(numUnits,numConds);
		
		//Initialize sample totals
		for(int j=0; j<numSamples; j++){
			double total =0;
			for(int i=0; i<numUnits; i++)
				total+= counts.get(i,j);
			totals[j]=total;
			scaling[j]=1;
		}
		//Initialize the unit -> index map
		unitToIndex = new HashMap<String, Integer>();
		for(int i=0; i<numUnits; i++){
			String un = unitNames[i];
			unitToIndex.put(un, i);
		}
	}
	
	/**
	 * Constructor: initialize data from events
	 * @param expts
	 */
	public CountsDataset(ExperimentManager manager, List<BindingEvent> events, int focalCond){
		this.focalCondition = focalCond;
		this.numConds = manager.getNumConditions();
		//Count points & samples
		int sampleCount = manager.getExperimentSet().getReplicates().size();
		int numPoints=events.size();
		
		//Set up counts array
        this.counts = new Matrix(numPoints,sampleCount);
        this.unitNames = new String[numPoints];
        int d=0;
        for(BindingEvent be : events){
			for(ExperimentCondition c : manager.getExperimentSet().getConditions()){
				for(ControlledExperiment rep : c.getReplicates()){
					counts.set(d, rep.getIndex(), be.getRepSigHits(rep));
				}
			}
			unitNames[d] = be.getPoint().getLocationString();
			d++;
        }
        
        //Name translators & design
        this.condToName = new HashMap<Integer, String>();
        this.sampleToExptName = new HashMap<Integer, Pair<String,String>>();
        this.design = new int[sampleCount];
        for(ExperimentCondition c : manager.getExperimentSet().getConditions()){
        	for(ControlledExperiment rep : c.getReplicates()){
        		design[rep.getIndex()] = c.getIndex();
        		sampleToExptName.put(rep.getIndex(), new Pair<String,String>(c.getName(), rep.getRepName()));
        	}condToName.put(c.getIndex(), c.getName());
        }
        
		this.numSamples = counts.getColumnDimension();
		this.numUnits = counts.getRowDimension();
		this.totals = new double[numSamples];
		this.scaling = new double[numSamples];
		this.DEpval = new Matrix(numUnits,numConds);
		this.condFold = new Matrix(numUnits,numConds);
		this.condMean = new Matrix(numUnits,numConds);
		this.condRawVar = new Matrix(numUnits,numConds);
		
		//Initialize sample totals
		for(int j=0; j<numSamples; j++){
			double total =0;
			for(int i=0; i<numUnits; i++)
				total+= counts.get(i,j);
			totals[j]=total;
			scaling[j]=1;
		}
		//Initialize the unit -> index map
		unitToIndex = new HashMap<String, Integer>();
		for(int i=0; i<numUnits; i++){
			String un = unitNames[i];
			unitToIndex.put(un, i);
		}
	}
	
	/**
	 * Calculcate the scaled mean and fold for all points in reference to the focal condition
	 */
	public void calcScMeanAndFold(){
		
		for(int i=0; i<numUnits; i++){
			//Average for focal condition
			double focalScTotal=0;
			double focalSamps=0;
			for(int s=0; s<numSamples; s++){
				if(design[s]==focalCondition){
					focalScTotal+= counts.get(i,s)/scaling[s];
					focalSamps++;
				}
			}
			double focalScMean = focalScTotal/focalSamps;
			condMean.set(i, focalCondition, focalScMean);
			condFold.set(i, focalCondition, 1.0);
			//Averages for other conditions
			for(int c=0; c<numConds; c++){
				if(c!=focalCondition){
					double otherScTotal=0;
					double otherSamps=0;
					for(int s=0; s<numSamples; s++){
						if(design[s]==c){
							otherScTotal+= counts.get(i,s)/scaling[s];
							otherSamps++;
						}
					}
					double otherScMean = otherScTotal/otherSamps;
					condMean.set(i, c, (otherScMean+focalScMean)/2);
					double fold = focalScMean>0 ? (otherScMean>0 ? (otherScMean/focalScMean) : 1/focalScMean) : otherScMean;
					condFold.set(i, c, fold);
				}
			}
		}
	}
	
	/**
	 * Update the provided events with differential expression values in the current dataset.  
	 * @param events
	 * @return
	 */
	public List<BindingEvent> updateEvents( List<BindingEvent> events, ExperimentManager manager){
		
		for(int e=0; e<events.size(); e++){
			BindingEvent ev = events.get(e);
			
			ExperimentCondition ref = manager.getExperimentSet().getIndexedCondition(focalCondition); 
			for(ExperimentCondition c : manager.getExperimentSet().getConditions()){
				if(c!=ref){
					//Only update the p-value since this should only be called after EdgeR as EdgeR's fold and mean are weird.
					ev.setInterCondP(c, ref, DEpval.get(e, c.getIndex()));
					ev.setInterCondFold(c, ref, condFold.get(e, c.getIndex()));
					ev.setInterCondScMean(c, ref, condMean.get(e, c.getIndex()));
				}
			}
			
		}
		return events;
	}
	
	/**
	 * Get the counts matrix
	 * @return
	 */
	public Matrix getCounts(){return counts;}
	
	/**
	 * Get an individual count
	 * @param u Unit index (i.e. gene or peak)
	 * @param s Sample index
	 * @return
	 */
	public double getCount(int u, int s){return counts.get(u, s);}
	
	/**
	 * Get the sum of counts for all samples
	 * @return
	 */
	public double[] getTotals(){return totals;}
	
	/**
	 * Get the sum of counts for sample s
	 * @param s
	 * @return
	 */
	public double getTotal(int s){return totals[s];}
	
	/**
	 * Get the name of unit u
	 * @param u
	 * @return
	 */
	public String getUnitName(int u){return unitNames[u];}
	
	/**
	 * Get the name of an experiment (Condition:Replicate)
	 * @param i
	 * @return
	 */
	public Pair<String,String> getExptName(int i){ return sampleToExptName.get(i);}
	
	/**
	 * Get the name of a condition (Condition:Replicate)
	 * @param i
	 * @return
	 */
	public String getCondName(int i){ return condToName.get(i);}
	
	/**
	 * Set the scaling factors array
	 * @param s double[]
	 */
	public void setScalingFactors(double[] s){
		if(s.length==scaling.length)
			scaling = s;
		else{
			System.err.println("Error: scaling array has incompatible length."); System.exit(1);
		}
	}
	/**
	 * Get the scaling array for all samples
	 * @return
	 */
	public double[] getScalingFactors(){return scaling;}
	/**
	 * Get the scaling factor for sample s
	 * @param s
	 * @return
	 */
	public double getScalingFactor(int s){return scaling[s];}
	
	/**
	 * Set the index of the condition that all comparisons are in reference to.
	 * @param x
	 */
	public void setFocalCondition(int x){focalCondition=x;}
	
	/**
	 * Get the focal condition
	 * @param x
	 * @return
	 */
	public int getFocalCondition(){return focalCondition;}
	
	/**
	 * Get the number of samples
	 * @return
	 */
	public int getNumSamples(){return numSamples;}
	
	/**
	 * Get the number of units (i.e. genes or peaks)
	 * @return
	 */
	public int getNumUnits(){return numUnits;}
	
	/**
	 * Get the number of conditions
	 * @return
	 */
	public int getNumConditions(){return numConds;}
	
	/**
	 * Get the design array (all samples)
	 * @return
	 */
	public int[] getDesignArray(){return design;}
	
	/**
	 * Get the condition corresponding to sample s
	 * @param i
	 * @return
	 */
	public int getSampleCondition(int s){return design[s];}
	
	
	/**
	 * Get the differential expression p-value matrix
	 * @return
	 */
	public Matrix getDEpvals(){return DEpval;}
	
	/**
	 * Get an individual p-value
	 * @param u Unit index (i.e. gene or peak)
	 * @param c Condition index
	 * @return
	 */
	public double getDEpval(int u, int c){return DEpval.get(u, c);}
	
	/**
	 * Set the DE p-values matrix
	 * @param p
	 */
	public void setDEpvals(Matrix p){DEpval=p;}
	
	/**
	 * Set an individual p-value
	 * @param u
	 * @param c
	 * @param p
	 */
	public void setDEpval(int u, int c, double p){DEpval.set(u, c,p);}
	
	
	/**
	 * Get the differential expression p-value matrix
	 * @return
	 */
	public Matrix getCondMeans(){return condMean;}
	
	/**
	 * Get an individual p-value
	 * @param u Unit index (i.e. gene or peak)
	 * @param c Condition index
	 * @return
	 */
	public double getCondMean(int u, int c){return condMean.get(u, c);}
	
	/**
	 * Get the differential expression fold matrix
	 * @return
	 */
	public Matrix getCondFolds(){return condFold;}
	
	/**
	 * Get an individual fold value
	 * @param u Unit index (i.e. gene or peak)
	 * @param c Condition index
	 * @return
	 */
	public double getCondFold(int u, int c){return condFold.get(u, c);}
	
	/**
	 * Set the conditions mean matrix
	 * @param c
	 */
	public void setCondMeans(Matrix c){condMean=c;}
	
	/**
	 * Set an individual p-value
	 * @param u
	 * @param c
	 * @param p
	 */
	public void setCondMean(int u, int c, double p){condMean.set(u, c,p);}
	
	/**
	 * Set the conditions mean matrix
	 * @param c
	 */
	public void setCondFolds(Matrix c){condFold=c;}
	
	/**
	 * Set an individual p-value
	 * @param u
	 * @param c
	 * @param p
	 */
	public void setCondFold(int u, int c, double f){condFold.set(u, c, f);}
	
	/**
	 * Get the differential expression p-value matrix
	 * @return
	 */
	public Matrix getCondRawVar(){return condRawVar;}
	
	/**
	 * Get an individual p-value
	 * @param u Unit index (i.e. gene or peak)
	 * @param c Condition index
	 * @return
	 */
	public double getCondRawVar(int u, int c){return condRawVar.get(u, c);}
	
	/**
	 * Set the conditions vars matrix
	 * @param c
	 */
	public void setCondRawVars(Matrix c){condRawVar=c;}
	
	/**
	 * Set an individual p-value
	 * @param u
	 * @param c
	 * @param p
	 */
	public void setCondRawVar(int u, int c, double p){condRawVar.set(u, c,p);}
	
	/**
	 * Return the scaling factors
	 * If tableFormat is true, print on a single line
	 * @param tableFormat
	 * @return
	 */
	public String scalingFactorString(boolean tableFormat){
		String str="";
		if(tableFormat){
			Pair<String,String> name = sampleToExptName.get(focalCondition);
			str = name.car()+":"+name.cdr();
			for(int s=0; s<numSamples; s++)
				str = str+"\t"+scaling[s];
		}else{
			for(int s=0; s<numSamples; s++){
				Pair<String,String> name = sampleToExptName.get(s);
				str = str+name.car()+":"+name.cdr()+"\t"+scaling[s]+"\n";
			}
		}
		return str;
	}
	
	/**
	 * Return the experiment names (condition:replicate)
	 * If tableFormat is true, print on a single line
	 * @param tableFormat
	 * @return
	 */
	public String experimentNamesString(boolean tableFormat){
		String str="";
		if(tableFormat){
			for(int s=0; s<numSamples; s++){
				Pair<String,String> name = sampleToExptName.get(s);
				str = str+"\t"+name.car()+":"+name.cdr();
			}
		}else{
			for(int s=0; s<numSamples; s++){
				Pair<String,String> name = sampleToExptName.get(s);
				str = str+name.car()+":"+name.cdr()+"\n";
			}
		}
		return str;
	}
	
	/**
	 * Return the index number for a given unit name
	 * @param s
	 * @return
	 */
	public Integer getUnitID(String s){
		if(unitToIndex.containsKey(s))
			return unitToIndex.get(s);
		else
			return -Integer.MAX_VALUE;
	}
	
	/**
	 * Scatters of each sample in the focal condition against one another 
	 * @param rasterImage
	 */
	public void savePairwiseFocalSampleMAPlots(String directory, boolean rasterImage){
		double A_min=1;
		int ref = 0;
		//Set one sample as the reference (the deepest sequenced sample in the focal condition)
		double maxTotal=0;
		for(int s=0; s<numSamples; s++)
			if(design[s] == focalCondition)
				if(totals[s]>maxTotal){
					ref=s; maxTotal=totals[s];
				}
		
		//Scale all focal condition samples against the reference
		for(int x=0; x<numSamples; x++){
			if(design[x]==focalCondition && x!=ref){
				List<Pair<Double,Double>> highlightMA = new ArrayList<Pair<Double,Double>>();
				List<Pair<Double,Double>> otherMA = new ArrayList<Pair<Double,Double>>();
				for(int d=0; d<numUnits; d++){
					double fold=0, avg=0;
					if(counts.get(d,x)>0 || counts.get(d,ref)>0)
						avg = ((counts.get(d,x)/scaling[x])+(counts.get(d,ref)/scaling[ref]))/2;
					if(counts.get(d,x)>0 && counts.get(d,ref)>0)
						fold = Math.log((counts.get(d,x)/scaling[x])/(counts.get(d,ref)/scaling[ref]))/LOG_2;
					else if(counts.get(d,x)>0)
						fold = Math.log((counts.get(d,x)/scaling[x]))/LOG_2;
					else if(counts.get(d,ref)>0)
						fold = Math.log((counts.get(d,ref)/scaling[ref]))/LOG_2;
					if(avg<A_min)
						avg = A_min;
					otherMA.add(new Pair<Double,Double>(fold,avg));
				}
					
				//Make the MA matrices
				Matrix maMatrixHighlight = new Matrix(highlightMA.size(),2);
				Matrix maMatrixOther = new Matrix(otherMA.size(),2);
				int count=0;
				for(Pair<Double,Double> v : highlightMA){
					maMatrixHighlight.set(count, 0, v.cdr());
					maMatrixHighlight.set(count, 1, v.car());
					count++;
				}
				count=0;
				for(Pair<Double,Double> v : otherMA){
					maMatrixOther.set(count, 0, v.cdr());
					maMatrixOther.set(count, 1, v.car());
					count++;
				}
				
				//Image name
				Pair<String,String> refName = getExptName(ref);
				Pair<String,String> currName = getExptName(x);
				String fileName = currName.car()+"-"+currName.cdr()+"_vs_"+refName.car()+"-"+refName.cdr()+".MA";
				if(rasterImage)
					fileName = directory+fileName+".png";
				else
					fileName = directory+fileName+".svg";
				//Generate image
				ScatterPlotMaker plotter = new ScatterPlotMaker(currName.car()+":"+currName.cdr()+" vs "+refName.car()+":"+refName.cdr()+" MA plot");
				plotter.saveMAplot(maMatrixOther, maMatrixHighlight, 0.0, fileName, rasterImage);
			}
		}
	}
	
	/**
	 * Scatters of each condition against the focal condition 
	 * @param rasterImage
	 */
	public void savePairwiseConditionMAPlots(double pValThreshold, String directory, boolean rasterImage){
		double A_min=1;
		
		//Scale all focal condition samples against the reference
		for(int c=0; c<numConds; c++){
			if(c!=focalCondition){
				List<Pair<Double,Double>> highlightMA = new ArrayList<Pair<Double,Double>>();
				List<Pair<Double,Double>> otherMA = new ArrayList<Pair<Double,Double>>();
				
				for(int d=0; d<numUnits; d++){
					double fold=condFold.get(d,c), avg = condMean.get(d,c);
					if(avg<A_min)
						avg = A_min;
					
					if(getDEpval(d, c)<pValThreshold)
						highlightMA.add(new Pair<Double,Double>(fold,avg));
					else
						otherMA.add(new Pair<Double,Double>(fold,avg));
				}
					
				//Make the MA matrices
				Matrix maMatrixHighlight = new Matrix(highlightMA.size(),2);
				Matrix maMatrixOther = new Matrix(otherMA.size(),2);
				int count=0;
				for(Pair<Double,Double> v : highlightMA){
					maMatrixHighlight.set(count, 0, v.cdr());
					maMatrixHighlight.set(count, 1, v.car());
					count++;
				}
				count=0;
				for(Pair<Double,Double> v : otherMA){
					maMatrixOther.set(count, 0, v.cdr());
					maMatrixOther.set(count, 1, v.car());
					count++;
				}
				
				//Image name
				String refName = getCondName(focalCondition);
				String currName = getCondName(c);
				String fileName = currName+"_vs_"+refName+".MA";
				if(rasterImage)
					fileName = directory+fileName+".png";
				else
					fileName = directory+fileName+".svg";
				//Generate image
				ScatterPlotMaker plotter = new ScatterPlotMaker(currName+" vs "+refName+" MA plot");
				plotter.saveMAplot(maMatrixOther, maMatrixHighlight, 0.0, fileName, rasterImage);
			}
		}
	}
	
	/**
	 * XY scatters of each sample in the focal condition against one another 
	 * @param rasterImage
	 */
	public void savePairwiseFocalSampleXYPlots(String directory, boolean rasterImage){
		int ref = 0;
		//Set one sample as the reference (the deepest sequenced sample in the focal condition)
		double maxTotal=0;
		for(int s=0; s<numSamples; s++)
			if(design[s] == focalCondition)
				if(totals[s]>maxTotal){
					ref=s; maxTotal=totals[s];
				}
		
		//Scale all focal condition samples against the reference
		for(int x=0; x<numSamples; x++){
			if(design[x]==focalCondition && x!=ref){
				List<Pair<Double,Double>> highlightMA = new ArrayList<Pair<Double,Double>>();
				List<Pair<Double,Double>> otherMA = new ArrayList<Pair<Double,Double>>();
				for(int d=0; d<numUnits; d++){
					otherMA.add(new Pair<Double,Double>(counts.get(d,x),counts.get(d,ref)));
				}
					
				//Make the MA matrices
				Matrix maMatrixHighlight = new Matrix(highlightMA.size(),2);
				Matrix maMatrixOther = new Matrix(otherMA.size(),2);
				int count=0;
				for(Pair<Double,Double> v : highlightMA){
					maMatrixHighlight.set(count, 0, v.car());
					maMatrixHighlight.set(count, 1, v.cdr());
					count++;
				}
				count=0;
				for(Pair<Double,Double> v : otherMA){
					maMatrixOther.set(count, 0, v.car());
					maMatrixOther.set(count, 1, v.cdr());
					count++;
				}
				
				//Image name
				Pair<String,String> refName = getExptName(ref);
				Pair<String,String> currName = getExptName(x);
				String fileName = currName.car()+"-"+currName.cdr()+"_vs_"+refName.car()+"-"+refName.cdr()+".XY";
				if(rasterImage)
					fileName = directory+fileName+".png";
				else
					fileName = directory+fileName+".svg";
				//Generate image
				ScatterPlotMaker plotter = new ScatterPlotMaker(currName.car()+":"+currName.cdr()+" vs "+refName.car()+":"+refName.cdr()+" XY plot");
				plotter.saveXYplot(maMatrixOther, maMatrixHighlight, currName.car()+"-"+currName.cdr(), refName.car()+"-"+refName.cdr(), fileName, rasterImage);
			}
		}
	}
	
	/**
	 * XY scatters of each condition against the focal condition 
	 * @param rasterImage
	 */
	public void savePairwiseConditionXYPlots(ExperimentManager man, double pValThreshold, String directory, boolean rasterImage){
		List<BindingEvent> events = man.getEvents();
		ExperimentCondition ref = man.getExperimentSet().getIndexedCondition(focalCondition); 
		for(ExperimentCondition c : man.getExperimentSet().getConditions()){
			if(c.getIndex()!=focalCondition){
				List<Pair<Double,Double>> highlightMA = new ArrayList<Pair<Double,Double>>();
				List<Pair<Double,Double>> otherMA = new ArrayList<Pair<Double,Double>>();
		
				for(int e=0; e<events.size(); e++){
					BindingEvent ev = events.get(e);			
					if(ev.getInterCondP(c, ref) < pValThreshold)
						highlightMA.add(new Pair<Double,Double>(ev.getCondSigHits(c), ev.getCondSigHits(ref)));
					else
						otherMA.add(new Pair<Double,Double>(ev.getCondSigHits(c), ev.getCondSigHits(ref)));
				}
				//Make the MA matrices
				Matrix maMatrixHighlight = new Matrix(highlightMA.size(),2);
				Matrix maMatrixOther = new Matrix(otherMA.size(),2);
				int count=0;
				for(Pair<Double,Double> v : highlightMA){
					maMatrixHighlight.set(count, 0, v.car());
					maMatrixHighlight.set(count, 1, v.cdr());
					count++;
				}
				count=0;
				for(Pair<Double,Double> v : otherMA){
					maMatrixOther.set(count, 0, v.car());
					maMatrixOther.set(count, 1, v.cdr());
					count++;
				}
				
				//Image name
				String refName = ref.getName();
				String currName = c.getName();
				String fileName = currName+"_vs_"+refName+".XY";
				if(rasterImage)
					fileName = directory+fileName+".png";
				else
					fileName = directory+fileName+".svg";
				//Generate image
				ScatterPlotMaker plotter = new ScatterPlotMaker(currName+" vs "+refName+" XY plot");
				plotter.saveXYplot(maMatrixOther, maMatrixHighlight, currName, refName, fileName, rasterImage);
			}
		}
	}
}
