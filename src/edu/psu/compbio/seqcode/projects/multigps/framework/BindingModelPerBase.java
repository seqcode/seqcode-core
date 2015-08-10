package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
	
	/**
	 * Generated using MetaMaker pile-ups of human TFIIB permanganate ChIP-seq (Pugh39735hg19)
	 * over initial binding event calls generated by Will Lai. 
	 */
	@SuppressWarnings("unchecked")
	public static final List<Pair<Integer, Double[]>> defaultPermChipSeqEmpiricalDistribution =
		Arrays.asList(
			 new Pair<Integer,Double[]>(-249, new Double[]{0.001884785, 0.001805877, 0.001697431, 0.002222478}),
			 new Pair<Integer,Double[]>(-247, new Double[]{0.002183365, 0.001993617, 0.002105556, 0.002034221}),
			 new Pair<Integer,Double[]>(-245, new Double[]{0.00210872, 0.002172417, 0.002068454, 0.002144037}),
			 new Pair<Integer,Double[]>(-243, new Double[]{0.002276672, 0.001984677, 0.001984974, 0.001919175}),
			 new Pair<Integer,Double[]>(-241, new Double[]{0.002183365, 0.002136657, 0.002448752, 0.002175414}),
			 new Pair<Integer,Double[]>(-239, new Double[]{0.002015414, 0.002020437, 0.002068454, 0.001777982}),
			 new Pair<Integer,Double[]>(-237, new Double[]{0.002220688, 0.001841637, 0.002133383, 0.002185872}),
			 new Pair<Integer,Double[]>(-235, new Double[]{0.002220688, 0.002324396, 0.00192932, 0.002112661}),
			 new Pair<Integer,Double[]>(-233, new Double[]{0.002276672, 0.002154537, 0.001855115, 0.002295689}),
			 new Pair<Integer,Double[]>(-231, new Double[]{0.002164704, 0.002002557, 0.002216863, 0.002112661}),
			 new Pair<Integer,Double[]>(-229, new Double[]{0.002743203, 0.002118777, 0.002374548, 0.002065597}),
			 new Pair<Integer,Double[]>(-227, new Double[]{0.002799186, 0.002109837, 0.00224469, 0.002212019}),
			 new Pair<Integer,Double[]>(-225, new Double[]{0.002575251, 0.002226056, 0.002133383, 0.002102203}),
			 new Pair<Integer,Double[]>(-223, new Double[]{0.002146043, 0.002297576, 0.002049903, 0.00228523}),
			 new Pair<Integer,Double[]>(-221, new Double[]{0.0024073, 0.002243936, 0.002550784, 0.00204468}),
			 new Pair<Integer,Double[]>(-219, new Double[]{0.002369978, 0.002234996, 0.00216121, 0.001997615}),
			 new Pair<Integer,Double[]>(-217, new Double[]{0.001884785, 0.002333336, 0.002133383, 0.002253854}),
			 new Pair<Integer,Double[]>(-215, new Double[]{0.001978092, 0.002190297, 0.002133383, 0.002070826}),
			 new Pair<Integer,Double[]>(-213, new Double[]{0.002537929, 0.002324396, 0.002318894, 0.002154496}),
			 new Pair<Integer,Double[]>(-211, new Double[]{0.00210872, 0.002386976, 0.001947871, 0.002212019}),
			 new Pair<Integer,Double[]>(-209, new Double[]{0.002351317, 0.002583656, 0.002309619, 0.002374129}),
			 new Pair<Integer,Double[]>(-207, new Double[]{0.002575251, 0.002136657, 0.001743809, 0.002248625}),
			 new Pair<Integer,Double[]>(-205, new Double[]{0.001978092, 0.002163477, 0.0020128, 0.001997615}),
			 new Pair<Integer,Double[]>(-203, new Double[]{0.002220688, 0.002547896, 0.002253965, 0.00203945}),
			 new Pair<Integer,Double[]>(-201, new Double[]{0.001940769, 0.002145597, 0.002151934, 0.002274771}),
			 new Pair<Integer,Double[]>(-199, new Double[]{0.001884785, 0.002243936, 0.002077729, 0.002159726}),
			 new Pair<Integer,Double[]>(-197, new Double[]{0.0024073, 0.002386976, 0.001892218, 0.002264313}),
			 new Pair<Integer,Double[]>(-195, new Double[]{0.002146043, 0.002234996, 0.00224469, 0.002133579}),
			 new Pair<Integer,Double[]>(-193, new Double[]{0.002295333, 0.002243936, 0.001938596, 0.002227707}),
			 new Pair<Integer,Double[]>(-191, new Double[]{0.001642189, 0.002226056, 0.002253965, 0.002055138}),
			 new Pair<Integer,Double[]>(-189, new Double[]{0.002052737, 0.002038317, 0.002383823, 0.001950551}),
			 new Pair<Integer,Double[]>(-187, new Double[]{0.002351317, 0.002449556, 0.002318894, 0.002175414}),
			 new Pair<Integer,Double[]>(-185, new Double[]{0.00255659, 0.002521076, 0.002374548, 0.002541469}),
			 new Pair<Integer,Double[]>(-183, new Double[]{0.001791479, 0.002485316, 0.002458028, 0.002102203}),
			 new Pair<Integer,Double[]>(-181, new Double[]{0.002127382, 0.002243936, 0.00216121, 0.002180643}),
			 new Pair<Integer,Double[]>(-179, new Double[]{0.001623528, 0.002306516, 0.002513681, 0.002295689}),
			 new Pair<Integer,Double[]>(-177, new Double[]{0.002369978, 0.002342276, 0.0020128, 0.002154496}),
			 new Pair<Integer,Double[]>(-175, new Double[]{0.002071398, 0.002449556, 0.002263241, 0.002384588}),
			 new Pair<Integer,Double[]>(-173, new Double[]{0.001791479, 0.002422736, 0.002105556, 0.002646056}),
			 new Pair<Integer,Double[]>(-171, new Double[]{0.00195943, 0.002369096, 0.001938596, 0.002253854}),
			 new Pair<Integer,Double[]>(-169, new Double[]{0.001436916, 0.002199237, 0.00241165, 0.00187734}),
			 new Pair<Integer,Double[]>(-167, new Double[]{0.002295333, 0.002243936, 0.002142658, 0.001866881}),
			 new Pair<Integer,Double[]>(-165, new Double[]{0.002519268, 0.002154537, 0.002068454, 0.00245257}),
			 new Pair<Integer,Double[]>(-163, new Double[]{0.002052737, 0.002243936, 0.002253965, 0.002295689}),
			 new Pair<Integer,Double[]>(-161, new Double[]{0.001940769, 0.002306516, 0.002522957, 0.002180643}),
			 new Pair<Integer,Double[]>(-159, new Double[]{0.002239349, 0.002324396, 0.002281792, 0.002170184}),
			 new Pair<Integer,Double[]>(-157, new Double[]{0.001996753, 0.002217116, 0.001920045, 0.002494405}),
			 new Pair<Integer,Double[]>(-155, new Double[]{0.002127382, 0.002208177, 0.002207587, 0.002295689}),
			 new Pair<Integer,Double[]>(-153, new Double[]{0.001940769, 0.002127717, 0.00209628, 0.002196331}),
			 new Pair<Integer,Double[]>(-151, new Double[]{0.002015414, 0.001823757, 0.00167888, 0.002212019}),
			 new Pair<Integer,Double[]>(-149, new Double[]{0.002649896, 0.002118777, 0.002355997, 0.00244734}),
			 new Pair<Integer,Double[]>(-147, new Double[]{0.002295333, 0.002190297, 0.001975698, 0.001950551}),
			 new Pair<Integer,Double[]>(-145, new Double[]{0.001978092, 0.002369096, 0.002189036, 0.002138808}),
			 new Pair<Integer,Double[]>(-143, new Double[]{0.001810141, 0.002646236, 0.00209628, 0.002055138}),
			 new Pair<Integer,Double[]>(-141, new Double[]{0.002220688, 0.001796937, 0.00241165, 0.00220156}),
			 new Pair<Integer,Double[]>(-139, new Double[]{0.002034075, 0.002315456, 0.001808738, 0.002253854}),
			 new Pair<Integer,Double[]>(-137, new Double[]{0.002351317, 0.002011497, 0.002291068, 0.002280001}),
			 new Pair<Integer,Double[]>(-135, new Double[]{0.001679512, 0.002127717, 0.002383823, 0.001981927}),
			 new Pair<Integer,Double[]>(-133, new Double[]{0.001791479, 0.002154537, 0.002374548, 0.00203945}),
			 new Pair<Integer,Double[]>(-131, new Double[]{0.001847463, 0.002011497, 0.00241165, 0.002133579}),
			 new Pair<Integer,Double[]>(-129, new Double[]{0.002034075, 0.002127717, 0.002346721, 0.002337524}),
			 new Pair<Integer,Double[]>(-127, new Double[]{0.001866124, 0.002011497, 0.00176236, 0.002274771}),
			 new Pair<Integer,Double[]>(-125, new Double[]{0.002071398, 0.002163477, 0.001957147, 0.001856423}),
			 new Pair<Integer,Double[]>(-123, new Double[]{0.001847463, 0.002190297, 0.002216863, 0.001981927}),
			 new Pair<Integer,Double[]>(-121, new Double[]{0.002052737, 0.002306516, 0.002207587, 0.002170184}),
			 new Pair<Integer,Double[]>(-119, new Double[]{0.002146043, 0.002226056, 0.002142658, 0.001981927}),
			 new Pair<Integer,Double[]>(-117, new Double[]{0.001492899, 0.002056197, 0.001938596, 0.001840735}),
			 new Pair<Integer,Double[]>(-115, new Double[]{0.001642189, 0.002324396, 0.001892218, 0.001804129}),
			 new Pair<Integer,Double[]>(-113, new Double[]{0.001604867, 0.002199237, 0.002105556, 0.001762294}),
			 new Pair<Integer,Double[]>(-111, new Double[]{0.002295333, 0.002422736, 0.002374548, 0.001872111}),
			 new Pair<Integer,Double[]>(-109, new Double[]{0.002146043, 0.001877397, 0.001920045, 0.001783212}),
			 new Pair<Integer,Double[]>(-107, new Double[]{0.001735496, 0.002333336, 0.002189036, 0.00187734}),
			 new Pair<Integer,Double[]>(-105, new Double[]{0.001978092, 0.002190297, 0.00224469, 0.001788441}),
			 new Pair<Integer,Double[]>(-103, new Double[]{0.002202027, 0.002252876, 0.002253965, 0.001762294}),
			 new Pair<Integer,Double[]>(-101, new Double[]{0.001884785, 0.002181357, 0.002207587, 0.001767524}),
			 new Pair<Integer,Double[]>(-99, new Double[]{0.002127382, 0.002351216, 0.001975698, 0.001683854}),
			 new Pair<Integer,Double[]>(-97, new Double[]{0.002295333, 0.002199237, 0.002430201, 0.00163156}),
			 new Pair<Integer,Double[]>(-95, new Double[]{0.002332655, 0.002404856, 0.002346721, 0.00179367}),
			 new Pair<Integer,Double[]>(-93, new Double[]{0.002034075, 0.002270756, 0.002189036, 0.001730918}),
			 new Pair<Integer,Double[]>(-91, new Double[]{0.002015414, 0.002351216, 0.002485855, 0.001563579}),
			 new Pair<Integer,Double[]>(-89, new Double[]{0.002295333, 0.002476376, 0.00249513, 0.001652478}),
			 new Pair<Integer,Double[]>(-87, new Double[]{0.001716834, 0.002351216, 0.00216121, 0.001845964}),
			 new Pair<Integer,Double[]>(-85, new Double[]{0.002164704, 0.002109837, 0.002485855, 0.001673395}),
			 new Pair<Integer,Double[]>(-83, new Double[]{0.002202027, 0.002306516, 0.002773398, 0.001464221}),
			 new Pair<Integer,Double[]>(-81, new Double[]{0.002369978, 0.002583656, 0.0028105, 0.001746606}),
			 new Pair<Integer,Double[]>(-79, new Double[]{0.00225801, 0.002315456, 0.002476579, 0.001710001}),
			 new Pair<Integer,Double[]>(-77, new Double[]{0.002836509, 0.002735636, 0.003079492, 0.001683854}),
			 new Pair<Integer,Double[]>(-75, new Double[]{0.003135089, 0.002851855, 0.002764122, 0.001929634}),
			 new Pair<Integer,Double[]>(-73, new Double[]{0.002276672, 0.002923375, 0.002578611, 0.001772753}),
			 new Pair<Integer,Double[]>(-71, new Double[]{0.002836509, 0.002655176, 0.003227901, 0.001783212}),
			 new Pair<Integer,Double[]>(-69, new Double[]{0.002761864, 0.002699876, 0.002856878, 0.002217249}),
			 new Pair<Integer,Double[]>(-67, new Double[]{0.002220688, 0.003182635, 0.003561822, 0.001893028}),
			 new Pair<Integer,Double[]>(-65, new Double[]{0.002817848, 0.002950195, 0.002986736, 0.002008074}),
			 new Pair<Integer,Double[]>(-63, new Double[]{0.00330304, 0.003316735, 0.003645302, 0.002154496}),
			 new Pair<Integer,Double[]>(-61, new Double[]{0.00300446, 0.003263095, 0.003543271, 0.002269542}),
			 new Pair<Integer,Double[]>(-59, new Double[]{0.002761864, 0.003629634, 0.00377516, 0.002274771}),
			 new Pair<Integer,Double[]>(-57, new Double[]{0.003265717, 0.003593874, 0.003255728, 0.002311377}),
			 new Pair<Integer,Double[]>(-55, new Double[]{0.003638942, 0.003727974, 0.003821538, 0.002379359}),
			 new Pair<Integer,Double[]>(-53, new Double[]{0.003788232, 0.003540234, 0.003728782, 0.002426423}),
			 new Pair<Integer,Double[]>(-51, new Double[]{0.003732248, 0.003620694, 0.003969947, 0.00261468}),
			 new Pair<Integer,Double[]>(-49, new Double[]{0.003209734, 0.004452113, 0.003515444, 0.00269835}),
			 new Pair<Integer,Double[]>(-47, new Double[]{0.003825555, 0.003978294, 0.004081254, 0.002630368}),
			 new Pair<Integer,Double[]>(-45, new Double[]{0.00345233, 0.004595153, 0.004619238, 0.002782019}),
			 new Pair<Integer,Double[]>(-43, new Double[]{0.003844216, 0.004166033, 0.004183285, 0.002818625}),
			 new Pair<Integer,Double[]>(-41, new Double[]{0.003937522, 0.004371653, 0.004582135, 0.00294413}),
			 new Pair<Integer,Double[]>(-39, new Double[]{0.0039002, 0.004711373, 0.00513867, 0.003111469}),
			 new Pair<Integer,Double[]>(-37, new Double[]{0.004441376, 0.004961692, 0.004359521, 0.003226515}),
			 new Pair<Integer,Double[]>(-35, new Double[]{0.004609327, 0.004729253, 0.005240701, 0.003472295}),
			 new Pair<Integer,Double[]>(-33, new Double[]{0.004702633, 0.005337172, 0.005593173, 0.004298534}),
			 new Pair<Integer,Double[]>(-31, new Double[]{0.005243809, 0.005077912, 0.005083016, 0.004110277}),
			 new Pair<Integer,Double[]>(-29, new Double[]{0.005691679, 0.004505753, 0.004628513, 0.003864497}),
			 new Pair<Integer,Double[]>(-27, new Double[]{0.006046243, 0.004943812, 0.004841851, 0.004601837}),
			 new Pair<Integer,Double[]>(-25, new Double[]{0.005654356, 0.005569611, 0.005333457, 0.004638442}),
			 new Pair<Integer,Double[]>(-23, new Double[]{0.005579711, 0.005355052, 0.005639551, 0.004560002}),
			 new Pair<Integer,Double[]>(-21, new Double[]{0.005803646, 0.006007671, 0.006353771, 0.004554772}),
			 new Pair<Integer,Double[]>(-19, new Double[]{0.0063075, 0.00655301, 0.006288841, 0.004811011}),
			 new Pair<Integer,Double[]>(-17, new Double[]{0.006120888, 0.005971911, 0.007253502, 0.00539147}),
			 new Pair<Integer,Double[]>(-15, new Double[]{0.006400806, 0.007196689, 0.007105092, 0.006306608}),
			 new Pair<Integer,Double[]>(-13, new Double[]{0.007408513, 0.008564506, 0.0076338, 0.005684314}),
			 new Pair<Integer,Double[]>(-11, new Double[]{0.006232855, 0.010093244, 0.008691216, 0.005496057}),
			 new Pair<Integer,Double[]>(-9, new Double[]{0.006624741, 0.010298864, 0.00803265, 0.007713306}),
			 new Pair<Integer,Double[]>(-7, new Double[]{0.007725754, 0.010111124, 0.007912068, 0.007477984}),
			 new Pair<Integer,Double[]>(-5, new Double[]{0.008397559, 0.010674343, 0.008654114, 0.010150187}),
			 new Pair<Integer,Double[]>(-3, new Double[]{0.010860843, 0.009994904, 0.008505704, 0.018893677}),
			 new Pair<Integer,Double[]>(-1, new Double[]{0.010562263, 0.009458505, 0.006706242, 0.048115339}),
			 new Pair<Integer,Double[]>(1, new Double[]{0.009871797, 0.008189027, 0.006873203, 0.076751313}),
			 new Pair<Integer,Double[]>(3, new Double[]{0.010767537, 0.009798225, 0.007429738, 0.028970653}),
			 new Pair<Integer,Double[]>(5, new Double[]{0.009106686, 0.009467445, 0.007299879, 0.015050097}),
			 new Pair<Integer,Double[]>(7, new Double[]{0.008136302, 0.008296307, 0.007763658, 0.009371013}),
			 new Pair<Integer,Double[]>(9, new Double[]{0.007296546, 0.009029386, 0.006696967, 0.008612755}),
			 new Pair<Integer,Double[]>(11, new Double[]{0.006997966, 0.008546626, 0.006687691, 0.008063673}),
			 new Pair<Integer,Double[]>(13, new Double[]{0.006251516, 0.007607928, 0.006233188, 0.006766791}),
			 new Pair<Integer,Double[]>(15, new Double[]{0.006232855, 0.006821209, 0.006149708, 0.007033489}),
			 new Pair<Integer,Double[]>(17, new Double[]{0.005878291, 0.00658877, 0.006233188, 0.00628569}),
			 new Pair<Integer,Double[]>(19, new Double[]{0.005169164, 0.005748411, 0.005852889, 0.00572092}),
			 new Pair<Integer,Double[]>(21, new Double[]{0.004758617, 0.005694771, 0.005296355, 0.005103855}),
			 new Pair<Integer,Double[]>(23, new Double[]{0.004534682, 0.005650071, 0.004851127, 0.004450185}),
			 new Pair<Integer,Double[]>(25, new Double[]{0.004553343, 0.004630913, 0.004943883, 0.004277616}),
			 new Pair<Integer,Double[]>(27, new Double[]{0.005169164, 0.005033212, 0.004637789, 0.003901102}),
			 new Pair<Integer,Double[]>(29, new Double[]{0.004105473, 0.004711373, 0.004313143, 0.003691928}),
			 new Pair<Integer,Double[]>(31, new Double[]{0.004348069, 0.004318013, 0.004443002, 0.003895873}),
			 new Pair<Integer,Double[]>(33, new Double[]{0.003508314, 0.004648793, 0.003932845, 0.003911561}),
			 new Pair<Integer,Double[]>(35, new Double[]{0.003844216, 0.004076634, 0.004007049, 0.0035089}),
			 new Pair<Integer,Double[]>(37, new Double[]{0.003974845, 0.004398473, 0.003645302, 0.003069634}),
			 new Pair<Integer,Double[]>(39, new Double[]{0.004292086, 0.003710094, 0.003357759, 0.00310624}),
			 new Pair<Integer,Double[]>(41, new Double[]{0.004236102, 0.003325675, 0.003348483, 0.002886607}),
			 new Pair<Integer,Double[]>(43, new Double[]{0.003340362, 0.003012775, 0.003496893, 0.002593762}),
			 new Pair<Integer,Double[]>(45, new Double[]{0.003470991, 0.003397195, 0.003070216, 0.002395047}),
			 new Pair<Integer,Double[]>(47, new Double[]{0.002817848, 0.003352495, 0.003237177, 0.002829084}),
			 new Pair<Integer,Double[]>(49, new Double[]{0.002575251, 0.003191575, 0.002912531, 0.002076056}),
			 new Pair<Integer,Double[]>(51, new Double[]{0.002836509, 0.002985955, 0.003033114, 0.002342753}),
			 new Pair<Integer,Double[]>(53, new Double[]{0.002369978, 0.002592596, 0.002532233, 0.002557157}),
			 new Pair<Integer,Double[]>(55, new Double[]{0.002743203, 0.003173695, 0.002634264, 0.002227707}),
			 new Pair<Integer,Double[]>(57, new Double[]{0.002034075, 0.002664116, 0.002170485, 0.002028992}),
			 new Pair<Integer,Double[]>(59, new Double[]{0.003041783, 0.002395916, 0.002235414, 0.00203945}),
			 new Pair<Integer,Double[]>(61, new Double[]{0.002313994, 0.001966797, 0.002207587, 0.002055138}),
			 new Pair<Integer,Double[]>(63, new Double[]{0.002313994, 0.002074077, 0.00209628, 0.001898258}),
			 new Pair<Integer,Double[]>(65, new Double[]{0.002444623, 0.002342276, 0.002235414, 0.001940092}),
			 new Pair<Integer,Double[]>(67, new Double[]{0.002164704, 0.001984677, 0.002040627, 0.001673395}),
			 new Pair<Integer,Double[]>(69, new Double[]{0.002313994, 0.002252876, 0.001873667, 0.001678624}),
			 new Pair<Integer,Double[]>(71, new Double[]{0.001530222, 0.001931037, 0.002105556, 0.001626331}),
			 new Pair<Integer,Double[]>(73, new Double[]{0.001436916, 0.001850577, 0.001808738, 0.001479909}),
			 new Pair<Integer,Double[]>(75, new Double[]{0.002202027, 0.001537678, 0.001743809, 0.001532202}),
			 new Pair<Integer,Double[]>(77, new Double[]{0.001847463, 0.001877397, 0.002022076, 0.001417156}),
			 new Pair<Integer,Double[]>(79, new Double[]{0.001698173, 0.001689657, 0.001957147, 0.00139101}),
			 new Pair<Integer,Double[]>(81, new Double[]{0.001604867, 0.001600257, 0.001558297, 0.001255046}),
			 new Pair<Integer,Double[]>(83, new Double[]{0.001399593, 0.001707537, 0.001400612, 0.001328257}),
			 new Pair<Integer,Double[]>(85, new Double[]{0.001642189, 0.001555558, 0.001688155, 0.001207982}),
			 new Pair<Integer,Double[]>(87, new Double[]{0.001492899, 0.001439338, 0.001511919, 0.001286422}),
			 new Pair<Integer,Double[]>(89, new Double[]{0.001250303, 0.001698597, 0.001586124, 0.001187065}),
			 new Pair<Integer,Double[]>(91, new Double[]{0.001343609, 0.001251598, 0.001456266, 0.001072019}),
			 new Pair<Integer,Double[]>(93, new Double[]{0.001231642, 0.001412518, 0.00167888, 0.001134771}),
			 new Pair<Integer,Double[]>(95, new Double[]{0.001250303, 0.001582377, 0.001344959, 0.001192294}),
			 new Pair<Integer,Double[]>(97, new Double[]{0.001418254, 0.001501918, 0.001344959, 0.001150459}),
			 new Pair<Integer,Double[]>(99, new Double[]{0.001474238, 0.001314178, 0.001242927, 0.00114523}),
			 new Pair<Integer,Double[]>(101, new Double[]{0.001418254, 0.001430398, 0.001270754, 0.001098166}),
			 new Pair<Integer,Double[]>(103, new Double[]{0.001362271, 0.001251598, 0.001558297, 0.001051101}),
			 new Pair<Integer,Double[]>(105, new Double[]{0.001548883, 0.001296298, 0.001715982, 0.001312569}),
			 new Pair<Integer,Double[]>(107, new Double[]{0.001231642, 0.001332058, 0.001224376, 0.001087707}),
			 new Pair<Integer,Double[]>(109, new Double[]{0.001399593, 0.001206898, 0.001326408, 0.001092936}),
			 new Pair<Integer,Double[]>(111, new Double[]{0.001212981, 0.001081738, 0.001224376, 0.001108624}),
			 new Pair<Integer,Double[]>(113, new Double[]{0.001511561, 0.001224778, 0.001354234, 0.000956973}),
			 new Pair<Integer,Double[]>(115, new Double[]{0.001287626, 0.001001278, 0.001344959, 0.000967432}),
			 new Pair<Integer,Double[]>(117, new Double[]{0.001268964, 0.001269478, 0.001493368, 0.00098312}),
			 new Pair<Integer,Double[]>(119, new Double[]{0.001101013, 0.001501918, 0.001011038, 0.001040643}),
			 new Pair<Integer,Double[]>(121, new Double[]{0.001156997, 0.001385698, 0.00136351, 0.001082477}),
			 new Pair<Integer,Double[]>(123, new Double[]{0.001362271, 0.001224778, 0.001289305, 0.001014496}),
			 new Pair<Integer,Double[]>(125, new Double[]{0.000933062, 0.001278418, 0.001474817, 0.000993578}),
			 new Pair<Integer,Double[]>(127, new Double[]{0.001287626, 0.001233718, 0.001372785, 0.000946514}),
			 new Pair<Integer,Double[]>(129, new Double[]{0.001306287, 0.001090678, 0.001261479, 0.000868074}),
			 new Pair<Integer,Double[]>(131, new Double[]{0.001380932, 0.001197958, 0.00128003, 0.001051101}),
			 new Pair<Integer,Double[]>(133, new Double[]{0.001212981, 0.001314178, 0.001335683, 0.001051101}),
			 new Pair<Integer,Double[]>(135, new Double[]{0.001082352, 0.001269478, 0.001419163, 0.000794863}),
			 new Pair<Integer,Double[]>(137, new Double[]{0.001231642, 0.001054918, 0.001242927, 0.000763487}),
			 new Pair<Integer,Double[]>(139, new Double[]{0.001287626, 0.001108558, 0.001558297, 0.001207982}),
			 new Pair<Integer,Double[]>(141, new Double[]{0.001119675, 0.000858239, 0.001344959, 0.000800092}),
			 new Pair<Integer,Double[]>(143, new Double[]{0.001175658, 0.001108558, 0.001159447, 0.00098312}),
			 new Pair<Integer,Double[]>(145, new Double[]{0.001343609, 0.001001278, 0.001233652, 0.000925597}),
			 new Pair<Integer,Double[]>(147, new Double[]{0.001436916, 0.001090678, 0.001224376, 0.000946514}),
			 new Pair<Integer,Double[]>(149, new Double[]{0.001138336, 0.001197958, 0.00113162, 0.000705964}),
			 new Pair<Integer,Double[]>(151, new Double[]{0.001250303, 0.001215838, 0.00096466, 0.000852386}),
			 new Pair<Integer,Double[]>(153, new Double[]{0.001250303, 0.001135378, 0.001242927, 0.000847156}),
			 new Pair<Integer,Double[]>(155, new Double[]{0.001380932, 0.001081738, 0.001085243, 0.000873303}),
			 new Pair<Integer,Double[]>(157, new Double[]{0.001324948, 0.001260538, 0.001029589, 0.000841927}),
			 new Pair<Integer,Double[]>(159, new Double[]{0.001380932, 0.001162198, 0.001233652, 0.000920367}),
			 new Pair<Integer,Double[]>(161, new Double[]{0.000914401, 0.001045978, 0.001242927, 0.000852386}),
			 new Pair<Integer,Double[]>(163, new Double[]{0.001082352, 0.001189018, 0.00104814, 0.000888991}),
			 new Pair<Integer,Double[]>(165, new Double[]{0.001343609, 0.001385698, 0.001215101, 0.000988349}),
			 new Pair<Integer,Double[]>(167, new Double[]{0.001156997, 0.001287358, 0.001354234, 0.000888991}),
			 new Pair<Integer,Double[]>(169, new Double[]{0.001324948, 0.001001278, 0.001094518, 0.000951743}),
			 new Pair<Integer,Double[]>(171, new Double[]{0.001101013, 0.001099618, 0.001382061, 0.000862844}),
			 new Pair<Integer,Double[]>(173, new Double[]{0.001063691, 0.000965518, 0.001344959, 0.000847156}),
			 new Pair<Integer,Double[]>(175, new Double[]{0.000858417, 0.000920819, 0.001428439, 0.000883762}),
			 new Pair<Integer,Double[]>(177, new Double[]{0.001343609, 0.001090678, 0.001029589, 0.000810551}),
			 new Pair<Integer,Double[]>(179, new Double[]{0.001175658, 0.001153258, 0.00136351, 0.000909909}),
			 new Pair<Integer,Double[]>(181, new Double[]{0.001156997, 0.001054918, 0.001419163, 0.000998808}),
			 new Pair<Integer,Double[]>(183, new Double[]{0.00089574, 0.001019158, 0.001437714, 0.000836698}),
			 new Pair<Integer,Double[]>(185, new Double[]{0.001138336, 0.000920819, 0.001168723, 0.000883762}),
			 new Pair<Integer,Double[]>(187, new Double[]{0.001455577, 0.001609197, 0.001326408, 0.001098166}),
			 new Pair<Integer,Double[]>(189, new Double[]{0.001418254, 0.001206898, 0.001140896, 0.000831468}),
			 new Pair<Integer,Double[]>(191, new Double[]{0.001287626, 0.000992338, 0.001233652, 0.000779175}),
			 new Pair<Integer,Double[]>(193, new Double[]{0.001250303, 0.001081738, 0.001372785, 0.000857615}),
			 new Pair<Integer,Double[]>(195, new Double[]{0.001082352, 0.001251598, 0.001521195, 0.000878532}),
			 new Pair<Integer,Double[]>(197, new Double[]{0.001343609, 0.001081738, 0.001437714, 0.001098166}),
			 new Pair<Integer,Double[]>(199, new Double[]{0.001250303, 0.001108558, 0.001242927, 0.001030184}),
			 new Pair<Integer,Double[]>(201, new Double[]{0.001268964, 0.001001278, 0.001326408, 0.001160918}),
			 new Pair<Integer,Double[]>(203, new Double[]{0.001754157, 0.000911879, 0.000983211, 0.001045872}),
			 new Pair<Integer,Double[]>(205, new Double[]{0.001268964, 0.001627077, 0.001289305, 0.000878532}),
			 new Pair<Integer,Double[]>(207, new Double[]{0.001156997, 0.001197958, 0.001187274, 0.000925597}),
			 new Pair<Integer,Double[]>(209, new Double[]{0.001324948, 0.001153258, 0.001317132, 0.001035413}),
			 new Pair<Integer,Double[]>(211, new Double[]{0.001138336, 0.001117498, 0.001521195, 0.001056331}),
			 new Pair<Integer,Double[]>(213, new Double[]{0.001250303, 0.001072798, 0.00136351, 0.001077248}),
			 new Pair<Integer,Double[]>(215, new Double[]{0.001119675, 0.001180078, 0.001400612, 0.001129542}),
			 new Pair<Integer,Double[]>(217, new Double[]{0.001007707, 0.001189018, 0.001539746, 0.001113854}),
			 new Pair<Integer,Double[]>(219, new Double[]{0.00119432, 0.001287358, 0.001344959, 0.00114523}),
			 new Pair<Integer,Double[]>(221, new Double[]{0.001380932, 0.001251598, 0.001307856, 0.001030184}),
			 new Pair<Integer,Double[]>(223, new Double[]{0.001343609, 0.001287358, 0.001400612, 0.000878532}),
			 new Pair<Integer,Double[]>(225, new Double[]{0.000933062, 0.001081738, 0.001372785, 0.000988349}),
			 new Pair<Integer,Double[]>(227, new Double[]{0.001306287, 0.001189018, 0.00136351, 0.001255046}),
			 new Pair<Integer,Double[]>(229, new Double[]{0.001250303, 0.001206898, 0.001382061, 0.000930826}),
			 new Pair<Integer,Double[]>(231, new Double[]{0.000951723, 0.001144318, 0.001317132, 0.001024955}),
			 new Pair<Integer,Double[]>(233, new Double[]{0.000933062, 0.001189018, 0.001224376, 0.000857615}),
			 new Pair<Integer,Double[]>(235, new Double[]{0.001138336, 0.001340998, 0.001066691, 0.000941285}),
			 new Pair<Integer,Double[]>(237, new Double[]{0.001101013, 0.001457218, 0.001335683, 0.001092936}),
			 new Pair<Integer,Double[]>(239, new Double[]{0.001436916, 0.001323118, 0.001224376, 0.001072019}),
			 new Pair<Integer,Double[]>(241, new Double[]{0.000989046, 0.001349938, 0.001344959, 0.000857615}),
			 new Pair<Integer,Double[]>(243, new Double[]{0.001511561, 0.001045978, 0.001474817, 0.001113854}),
			 new Pair<Integer,Double[]>(245, new Double[]{0.001492899, 0.001457218, 0.00113162, 0.000878532}),
			 new Pair<Integer,Double[]>(247, new Double[]{0.000653143, 0.001189018, 0.00136351, 0.001176606}),
			 new Pair<Integer,Double[]>(249, new Double[]{0.000970385, 0.001171138, 0.001270754, 0.001035413}));
	
}
