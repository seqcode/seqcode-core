package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;

/**
 * TagDistribution defines a (probabilistic) model of sequenced tag occurrences around a binding event.
 * The probability is directional (i.e. stranded). The probability distributions can be asymmetric on Watson & Crick strands. 
 * Given a signed distance from the event, the TagDistribution should return 
 * a relative probability of seeing a read at that distance (given that the tag is on Watson or Crick).
 * 
 * @author shaunmahony
 *
 */
public class TagDistribution {
	
	protected final double LOG2 = Math.log(2);
	protected int winSize;
	protected int left, right; //relative left & right positions
	protected double[] watsonData, crickData; //Data lanscape should be based on (typically tag counts)
	protected double[] watsonProbs, crickProbs; //Probability landscapes
	protected double[] watsonLogProbs, crickLogProbs; //Log probabilities
	protected int watsonSummit, crickSummit;		// relative positions of highest probs
	protected int influenceRange; //95% probability range (over both strands)
	protected double bgProb, logBgProb;
	
	public TagDistribution(int size){
		init(-(winSize/2), winSize-(winSize/2));
	}
	
	
	//Accessors
	public int getLeft(){return left;}
	public int getRight(){return right;}
	public int getWatsonSummit(){return watsonSummit;}
	public int getCrickSummit(){return crickSummit;}
	public int getInfluenceRange(){return influenceRange;}
	public double[] getWatsonProbabilities(){return  watsonProbs.clone();}
	public double[] getCrickProbabilities(){return  crickProbs.clone();}
	
	/**
	 *  Look up the probability corresponding to a distance and strand relationship to the central position.
	 *  Distance should be defined as (tag position - central position)
	 */
	public double probability(int distance, boolean watsonStrand){
		if(distance<left || distance>right){
			return(bgProb);
		}else{
			return(watsonStrand ? 
					watsonProbs[distance-left] : crickProbs[distance-left]);
		}
	}
	/**
	 *  Look up the log probability corresponding to a distance and strand relationship to the central position.
	 *  Distance should be defined as (tag position - central position)
	 */
	public double logProbability(int distance, boolean watsonStrand){
		if(distance<left || distance>right){
		  return(logBgProb);
		}else{
			return(watsonStrand ? 
					watsonLogProbs[distance-left] : crickLogProbs[distance-left]);
		}	  
	}
	
	/**
	 * Load data & make probabilities
	 * Input lists are pairs of relative positions and tag counts.
	 * Assumes the input lists are sorted according to increasing position.
	 * If the crick list is null, the mirror image of the watson list is used.   
	 * @param watsonTagDist: list of paired ints and doubles
	 * @param crickTagDist: list of paired ints and doubles. If null, uses mirror of watsonTagDist
	 * @throws Exception 
	 */
	protected void loadData(List<Pair<Integer, Double>> watsonTagDist, List<Pair<Integer, Double>> crickTagDist) throws Exception{		
		//Find left, right values first
		for(Pair<Integer, Double> p : watsonTagDist){
			if(p.car()<left){ left=p.car();}
			if(p.car()>right){right=p.car();}
		}
		for(Pair<Integer, Double> p : crickTagDist){
			if(p.car()<left){ left=p.car();}
			if(p.car()>right){right=p.car();}
		}
		
		//Initialize arrays
		init(left, right);
		
		//Populate the data array (assumes sorted)
		//WATSON
		int last=left-1;
		for(Pair<Integer, Double> p : watsonTagDist){
			int index = p.car();
			double val = p.cdr();
			//if list is not properly sorted, throw exception
			if(index-last<0)
				throw new Exception("Incorrectly sorted binding read distribution data!"); 
			//if unevenly spaced, linearly interpolate between values
			if(index-last>1){
				double lastVal=dataVal(last, true);
				double step = (val-lastVal)/(double)(index-last);
				for(int i=1; i<(index-last); i++){
					watsonData[(last+i)-left]=lastVal+(step*(double)i);
				}
			}
			watsonData[index-left]=val;
			last = p.car();
		}
		//CRICK
		if(crickTagDist==null){
			crickData = watsonData.clone();
		}else{
			last=left-1;
			for(Pair<Integer, Double> p : crickTagDist){
				int index = p.car();
				double val = p.cdr();
				//if list is not properly sorted, throw exception
				if(index-last<0)
					throw new Exception("Incorrectly sorted binding read distribution data!"); 
				//if unevenly spaced, linearly interpolate between values
				if(index-last>1){
					double lastVal=dataVal(last, false);
					double step = (val-lastVal)/(double)(index-last);
					for(int i=1; i<(index-last); i++){
						crickData[(last+i)-left]=lastVal+(step*(double)i);
					}
				}
				crickData[index-left]=val;
				last = p.car();
			}
		}
		makeProbabilities();
	}
	
	//Set a probability landscape according to the data. 
	protected void makeProbabilities(){
		double totalW=0, totalC=0, minProb=Double.MAX_VALUE;
		for(int i=left; i<=right; i++){
			totalW+=dataVal(i, true);
			totalC+=dataVal(i, false);
		}
		for(int i=left; i<=right; i++){
			watsonProbs[i-left] = dataVal(i, true)/totalW;
			crickProbs[i-left] = dataVal(i, false)/totalC;
			watsonLogProbs[i-left] = Math.log(watsonProbs[i-left])/LOG2;
			crickLogProbs[i-left] = Math.log(crickProbs[i-left])/LOG2;
			if(watsonProbs[i-left]<minProb){minProb = watsonProbs[i-left];}
			if(crickProbs[i-left]<minProb){minProb = crickProbs[i-left];}
		}
		Pair<Double, TreeSet<Integer>> wSorted = StatUtil.findMax(watsonProbs);
		Pair<Double, TreeSet<Integer>> cSorted = StatUtil.findMax(crickProbs);
		watsonSummit = wSorted.cdr().first()+left;
		crickSummit = cSorted.cdr().first()+left;
		
		bgProb = minProb/1000;
		logBgProb = Math.log(bgProb)/LOG2;

		updateInfluenceRange();
	}
	
	//Initialize the data structures
	protected void init(int left, int right){
		this.left = left;
		this.right = right;
		winSize = right-left+1;
		watsonData = new double [winSize];
		crickData = new double [winSize];
		watsonProbs = new double [winSize];
		crickProbs = new double [winSize];
		watsonLogProbs = new double [winSize];
		crickLogProbs = new double [winSize];
		for(int w=0; w<winSize; w++){
			watsonData[w]=0; crickData[w]=0;
			watsonProbs[w]=0; crickProbs[w]=0;
			watsonLogProbs[w]=Double.NEGATIVE_INFINITY; crickLogProbs[w]=Double.NEGATIVE_INFINITY;
		}
	}
	
	//Return a pair of distances corresponding to the central probability interval provided
	protected Pair<Integer,Integer> probIntervalDistances(double prob){
		double ends=(1-prob)/2;
		double probSum=0;
		boolean firstFound=false, secondFound=false;
		int first=left, second=right;
		for(int i=left; i<=right; i++){
			probSum+=probability(i, true)+probability(i, false);
			if(!firstFound && probSum>ends){
				firstFound=true;
				first=i;
			}else if(!secondFound && probSum>(1-ends)){
				secondFound=true;
				second=i;
			}
		}
		Pair<Integer,Integer> intervalDists = new Pair<Integer, Integer>(first, second);
		return intervalDists;
	}
	
	//Update the influence range
	protected void updateInfluenceRange(){
		Pair<Integer,Integer> intervals = probIntervalDistances(0.95);
		int longest = Math.max(Math.abs(intervals.car()), Math.abs(intervals.cdr()));
		influenceRange = longest*2;
	}
	
	//Look up the data corresponding to a distance
	protected double dataVal(int distance, boolean watsonStrand){
		if(distance<left || distance>right){
			return(0.0);
		}else{
			return(watsonStrand ? 
					watsonData[distance-left] : crickData[distance-left]);
		}
	}
}
