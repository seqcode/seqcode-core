package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.utils.Pair;

public class BindingModelPerBase extends BindingModel{
	protected double[][] data_pb;
	protected double[][] probs_pb;
	protected double[][] logProbs_pb;	

	protected List<Pair<Integer, Double[]>> empiricalDistributionPerBase; //Per base empirical distrib
	protected Map<Character, Integer> baseToInt;
	
	public BindingModelPerBase(File f, int minDist, int maxDist){
		super();
		initBaseToInt();
		min=0; max=0;
		fileName = f.getName();
		try {
			empiricalDistribution = new LinkedList<Pair<Integer,Double>>();
			empiricalDistributionPerBase = new LinkedList<Pair<Integer,Double[]>>();
			BufferedReader reader = new BufferedReader(new FileReader(f));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=5){
	              Integer dist = new Integer(words[0]);
	              //make sure the current data point is within the specified range
	              if ((dist.intValue() >= minDist) && (dist.intValue() <= maxDist)) {
	                Double[] curr = new Double[4];
	                for(int w=1; w<=4; w++){
	                	Double val = new Double(words[w]); 
	                	curr[w-1] = val;
	                
	                	//empiricalDistribution sums over all 4 bases
	                	Pair<Integer,Double> p = new Pair<Integer,Double>(dist, val);
		                if (p.cdr().doubleValue()>=0)	// should be non-negative value
		                	empiricalDistribution.add(p);
		                else {
		                	System.err.println("\nRead distribution file contains negative probability(count) value!"); 
		                	System.exit(1);
		                }
	                }
	                empiricalDistributionPerBase.add(new Pair<Integer, Double[]>(dist, curr));
	              }
	            }
	        }
	        loadData(empiricalDistribution);
	        loadDataPB(empiricalDistributionPerBase);
			makeProbabilities();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public BindingModelPerBase(File f) {
		this(f, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	
	public BindingModelPerBase(List<Pair<Integer, Double[]>> bindingDist){
		initBaseToInt();
		min=0; max=0;
		empiricalDistributionPerBase=bindingDist;
		empiricalDistribution = new LinkedList<Pair<Integer,Double>>();
		for(Pair<Integer,Double[]> ps : bindingDist){
			Pair<Integer,Double> p = new Pair<Integer,Double>(ps.car(), ps.cdr()[0]+ps.cdr()[1]+ps.cdr()[2]+ps.cdr()[3]);
			empiricalDistribution.add(p);
		}
		loadData(empiricalDistribution);
		loadDataPB(empiricalDistributionPerBase);
		makeProbabilities();
	}
	
	//Accessors
	public double[][] getProbabilitiesPerBase(){	return  probs_pb.clone();}
	public double[][] getLogProbabilitiesPerBase() { return logProbs_pb.clone();}
	public String getFileName() {
		return fileName;
	}
	public void setFileName(String fileName) {
		this.fileName = fileName;
	}	
	//Simple lookup
	private void initBaseToInt(){
		baseToInt = new HashMap<Character, Integer>();
		baseToInt.put('A', 0); baseToInt.put('a', 0);
		baseToInt.put('C', 1); baseToInt.put('c', 1);
		baseToInt.put('G', 2); baseToInt.put('g', 2);
		baseToInt.put('T', 3); baseToInt.put('t', 3);
		baseToInt.put('N', 0); baseToInt.put('n', 0);
	}
	
	//Load data
	protected void loadDataPB(List<Pair<Integer, Double[]>> bindingDist){
		//Assumes the list is sorted//
		
		//Find max, min values first
		for(Pair<Integer, Double[]> p : bindingDist){
			if(p.car()<min)
				min=p.car();
			if(p.car()>max)
				max=p.car();
		}
		//Initialize arrays
		data_pb = new double[(max-min)+1][4];
		probs_pb = new double[(max-min)+1][4];
		logProbs_pb = new double[(max-min)+1][4];
		for(int i=0; i<=(max-min); i++){
			data[i]=0; probs[i]=0; logProbs[i] = Double.NEGATIVE_INFINITY;
		}
		
		//Populate the data array (assumes sorted)
		int last=min-1;
		for(Pair<Integer, Double[]> p : bindingDist){
			int index = p.car();
			Double[] vals = p.cdr();
			//if list is not properly sorted (need to make this into an exception)
			if(index-last<0){
				System.err.println("Incorrectly sorted binding read distribution data!"); 
				System.exit(1);
			}
			for(int v=0; v<4; v++){
				//if unevenly spaced, smooth linearly between values
				if(index-last>1){
					double[] lastVals=dataPBVals(last);
					double step = (vals[v]-lastVals[v])/(double)(index-last);
					for(int i=1; i<(index-last); i++){
						data_pb[(last+i)-min][v]=lastVals[v]+(step*(double)i);
					}
				}
				data_pb[index-min][v]=vals[v];
			}
			
			
			last = p.car();
		}
	}
	
	//Set a probability landscape according to the data. 
	protected void makeProbabilities(){
		super.makeProbabilities();//Call parent method to do everything for the combined empirical distribution
		double[] totalVals={0,0,0,0}, minProbs={Double.MAX_VALUE,Double.MAX_VALUE,Double.MAX_VALUE,Double.MAX_VALUE};
		for(int i=min; i<=max; i++){
			for(int v=0; v<4; v++)
				totalVals[v] += dataPBVals(i)[v];
		}
		for(int i=min; i<=max; i++){
			for(int v=0; v<4; v++){
				probs_pb[i-min][v] = dataPBVals(i)[v]/totalVals[v]; 
				logProbs_pb[i-min][v] = Math.log(probs_pb[i-min][v])/LOG2;
				if(probs_pb[i-min][v]<minProbs[v])
					minProbs[v] = probs_pb[i-min][v];
			}
		}
		
		// update empiricalDistribution with normalized probability
		List<Pair<Integer, Double[]>> newDist = new ArrayList<Pair<Integer, Double[]>> ();
		for(int i=min; i<=max; i++){
			Double[] newArr = {probability(i, 'A'), probability(i, 'C'), probability(i, 'G'), probability(i, 'T')};
			newDist.add(new Pair<Integer, Double[]>(i, newArr));
		}
		empiricalDistributionPerBase=newDist;
	}
	
		
	/**
	 * Look up the probability corresponding to a distance
	 * @param distance : Distance should be defined as (Read position - Peak position)
	 * @param base : base that defines which component of the read distribution to use
	 * @return
	 */
	public double probability(int distance, char base){
		int baseIndex = baseToInt.get(base);
		if(distance<min || distance>max || baseIndex<0 || baseIndex>3){
			return(bgProb);
		}else{
			return(probs_pb[distance-min][baseIndex]);
		}
	}
	
	public double logProbability(int distance, char base) {
		int baseIndex = baseToInt.get(base);
		if(distance<min || distance>max || baseIndex<0 || baseIndex>3){
		  return(logBgProb);
		}else{
		  return(logProbs_pb[distance-min][baseIndex]);
		}	  
	}
	
	//Look up the data corresponding to a distance
	public double[] dataPBVals(int distance){
		if(distance<min || distance>max){
			double[] zeroArray = {0.0,0.0, 0.0, 0.0};
			return(zeroArray);
		}else{
			return(data_pb[distance-min]);
		}
	}
		
	//Print probs to a file
	public void printToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(int i=min; i<=max; i++){
				fout.write(i+"\t"+probability(i, 'A')+"\t"+probability(i, 'C')+"\t"+probability(i, 'G')+"\t"+probability(i, 'T')+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	public List<Pair<Integer, Double[]>> getEmpiricalPerBaseDistribution() {
		List<Pair<Integer, Double[]>> newDist = new ArrayList<Pair<Integer, Double[]>> ();
		for (Pair<Integer, Double[]> p: empiricalDistributionPerBase)
			newDist.add(p);
		return newDist;
	}
	
}
