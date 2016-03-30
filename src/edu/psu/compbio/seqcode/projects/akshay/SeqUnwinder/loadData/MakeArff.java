package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.loadData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;

import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;
import weka.core.converters.ConverterUtils.DataSource;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.RepeatMaskedRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.loadData.MakeArffFromFasta.MapKeyComparator;

public class MakeArff {
	
	/** The design file that contains the relationship between subgroups and labels */
	protected String design;
	protected List<String> subGroupNames = new ArrayList<String>();
	protected HashMap<String,Double> subGroupWeights = new HashMap<String,Double>();
	protected List<Point> peaks = new ArrayList<Point>();
	protected List<Region> regions = new ArrayList<Region>();
	protected GenomeConfig gcon;
	protected SequenceGenerator<Region> seqgen = null;
	protected int win = 150;
	protected String outArffFileName="out.arff";
	protected int Kmin=4;
	protected int Kmax=6;
	
	@SuppressWarnings("unchecked")
	public MakeArff(GenomeConfig gc) {
		gcon = gc;
		seqgen = gc.getSequenceGenerator();
	}
	
	/**
	 * Use this settor method if you only have subgroups (simple multi-logistic regression)
	 * @param labs
	 * @param randSeqs
	 */
	public void setLabelsSimple(List<String> labs, List<Region> randRegs){
		subGroupNames.addAll(labs);
		// Now make the design file
		StringBuilder designBuilder = new StringBuilder();
		//First get the header
		designBuilder.append("#Index"+"\t"+
				"Name"+"\t"+
				"isSubGroup" +"\t"+
				"AssignedLabs" + "\t"+
				"AssignedSGs" + "\t"+
				"Layer"+"\n");
		// Now subgroups
		HashMap<String,Integer> subgroupNames = new HashMap<String,Integer>();
		for(String s : subGroupNames){
			if(!subgroupNames.containsKey(s))
				subgroupNames.put(s, 1);
		}
		int index = 0;
		for(String s : subgroupNames.keySet()){
			designBuilder.append(index);designBuilder.append("\t"); // subgroup id
			designBuilder.append(s+"\t"); // subgroup 
			designBuilder.append(1);designBuilder.append("\t"); // subgroup indicator
			designBuilder.append("-"+"\t"); // Assigned labels
			designBuilder.append("-"+"\t"); // Assigned subgroups
			designBuilder.append(0);designBuilder.append("\n"); // Layer
			index++;
		}

		// Now add the random label subgroup
		designBuilder.append(index);designBuilder.append("\t"); 
		designBuilder.append("RootRandom"+"\t");
		designBuilder.append(1);designBuilder.append("\t");
		designBuilder.append("-"+"\t");
		designBuilder.append("-"+"\t");
		designBuilder.append(0);designBuilder.append("\n");

		designBuilder.deleteCharAt(designBuilder.length()-1);
		design = designBuilder.toString();
		
		// Now add the random regions to peaks and "RootRandom" label to subGroupNames
		regions.addAll(randRegs);
		for(Region re : randRegs){
			peaks.add(re.getMidpoint());
			subGroupNames.add("RootRandom");
		}

		// Finally, find the weights for each subgroup
		for(String s : subGroupNames){
			if(!subGroupWeights.containsKey(s))
				subGroupWeights.put(s, 1.0);
			else
				subGroupWeights.put(s, subGroupWeights.get(s)+1);
		}
		// Now get the sorted keys
		MapKeyComparator<Double> comp_d = new MapKeyComparator<Double>(subGroupWeights);
		List<String> groupWeightsKeyset = comp_d.getKeyList();
		Collections.sort(groupWeightsKeyset, comp_d);

		for(int i=0; i<groupWeightsKeyset.size()-1; i++){
			subGroupWeights.put(groupWeightsKeyset.get(i), 
					subGroupWeights.get(groupWeightsKeyset.get(groupWeightsKeyset.size()-1))/subGroupWeights.get(groupWeightsKeyset.get(i)));
		}
		subGroupWeights.put(groupWeightsKeyset.get(groupWeightsKeyset.size()-1), 1.0);
	}

	
	//Setters
	/**
	 * This setter methods does four things :-
	 * 1) Stores subgroup names in "subGroupNames". The format is "label1&label2"
	 * 2) Makes the design file for SeqUnwinder
	 * 3) Appends random regions to peaks and regions and also subGroupNames
	 * 4) Finally, calculates the weights for each subgroup instances.
	 * @param labs
	 * @param randRegs
	 */
	@SuppressWarnings("unchecked")
	public void setLabels(List<String> labs, List<Region> randRegs){
		HashMap<String,Integer> subgroupNames = new HashMap<String,Integer>();
		HashMap<String,Integer> labelNames = new HashMap<String,Integer>();
		int indSG=0;
		for(String s : labs){
			String[] labels = s.split(";");
			// Now generate the sub group name
			StringBuilder subgroupSB = new StringBuilder();
			for(String sname : labels){
				subgroupSB.append(sname);subgroupSB.append("&");
			}
			//String subgroup = labels[0]+"&"+labels[1];
			String subgroup = subgroupSB.toString();
			if(!subgroupNames.containsKey(subgroup)){
				subgroupNames.put(subgroup, indSG);
				indSG++;
			}
			for(String sname : labels){
				if(!labelNames.containsKey(sname)){
					labelNames.put(sname, 1);
				}
			}
			
			subGroupNames.add(subgroup);
		}
		
		// Now remove a label that has only one subgroup to it
		HashMap<String,Integer> tmpLabMap = new HashMap<String,Integer>(); 
		for(String sname : subgroupNames.keySet()){
			for(String piece : sname.split("&")){
				if(tmpLabMap.containsKey(piece)){
					tmpLabMap.put(piece,tmpLabMap.get(piece) + 1);
				}else{
					tmpLabMap.put(piece, 1);
				}
			}
		}
		
		for(String sname : tmpLabMap.keySet()){
			if(tmpLabMap.get(sname) ==1){
				labelNames.remove(sname);
			}
		}
		tmpLabMap.clear();
		
		// Now, adjest the indexes of the labels
		int indLab=0;
		for(String s : labelNames.keySet()){
			labelNames.put(s, indLab+indSG+1);
			indLab++;
		}
		
		// Now make the design file
		StringBuilder designBuilder = new StringBuilder();
		//First get the header
		designBuilder.append("#Index"+"\t"+
							 "Name"+"\t"+
							 "isSubGroup" +"\t"+
							 "AssignedLabs" + "\t"+
							 "AssignedSGs" + "\t"+
							 "Layer"+"\n");
		// Now subgroups
		int index = 0;
		MapKeyComparator<Integer> comp = new MapKeyComparator<Integer>(subgroupNames);
		List<String> subGs = comp.getKeyList();
		Collections.sort(subGs, comp);
		for(String s : subGs){
			// Check if this subgroup has parents
			boolean isRoot = true;
			for(String labname: labelNames.keySet()){
				if(s.startsWith(labname+"&") || s.endsWith("&"+labname) || s.contains("&"+labname+"&"))
					isRoot = false;
			}
			if(labelNames.size() == 0)
				isRoot = true;
			//If this subgroup is a root add the root tag to all the names in "subGroupNames"
			if(isRoot)
				addRootTag(s);
			designBuilder.append(index);designBuilder.append("\t"); // subgroup id
			designBuilder.append(s+"\t"); // subgroup 
			designBuilder.append(1);designBuilder.append("\t"); // subgroup indicator
			if(!isRoot){
				String[] assignedLabs = s.split("&");
				for(String aS : assignedLabs){
					designBuilder.append(labelNames.get(aS)+",");// Assigned labels
				}
				designBuilder.deleteCharAt(designBuilder.length()-1);
				designBuilder.append("\t");
			}else{
				designBuilder.append("-"+"\t");
			}
			designBuilder.append("-"+"\t"); // Assigned subgroups
			designBuilder.append(0);designBuilder.append("\n"); // Layer
			index++;
		}
		
		// Now add the random label subgroup
		designBuilder.append(index);designBuilder.append("\t"); 
		designBuilder.append("RootRandom"+"\t");
		designBuilder.append(1);designBuilder.append("\t");
		designBuilder.append("-"+"\t");
		designBuilder.append("-"+"\t");
		designBuilder.append(0);designBuilder.append("\n");
		
		index++;
		// Now Labels
		comp = new MapKeyComparator<Integer>(labelNames);
		List<String> LABs = comp.getKeyList();
		Collections.sort(LABs, comp);
		for(String s : LABs){
			designBuilder.append(index);designBuilder.append("\t"); // label id
			designBuilder.append("Root"+s+"\t"); // label name
			designBuilder.append(0);designBuilder.append("\t"); // subgroup indicator
			designBuilder.append("-"+"\t"); // Assigned labels
			
			for(String subS : subgroupNames.keySet()){ 
				if(subS.startsWith(s+"&") || subS.endsWith("&"+s) || subS.contains("&"+s+"&")){
					designBuilder.append(subgroupNames.get(subS)+","); // Assigned subgroups
				}
			}
			designBuilder.deleteCharAt(designBuilder.length()-1);
			designBuilder.append("\t");
			designBuilder.append(1);designBuilder.append("\n");
			index++;
		}
		designBuilder.deleteCharAt(designBuilder.length()-1);
		design = designBuilder.toString();
		
		// Now add the random regions to peaks and "RootRandom" label to subGroupNames
		regions.addAll(randRegs);
		for(Region re : randRegs){
			peaks.add(re.getMidpoint());
			subGroupNames.add("RootRandom");
		}
		
		// Finally, find the weights for each subgroup
		for(String s : subGroupNames){
			if(!subGroupWeights.containsKey(s))
				subGroupWeights.put(s, 1.0);
			else
				subGroupWeights.put(s, subGroupWeights.get(s)+1);
		}
		// Now get the sorted keys
		MapKeyComparator<Double> comp_d = new MapKeyComparator<Double>(subGroupWeights);
		List<String> groupWeightsKeyset = comp_d.getKeyList();
		Collections.sort(groupWeightsKeyset, comp_d);
		
		for(int i=0; i<groupWeightsKeyset.size()-1; i++){
			subGroupWeights.put(groupWeightsKeyset.get(i), 
					subGroupWeights.get(groupWeightsKeyset.get(groupWeightsKeyset.size()-1))/subGroupWeights.get(groupWeightsKeyset.get(i)));
		}
		subGroupWeights.put(groupWeightsKeyset.get(groupWeightsKeyset.size()-1), 1.0);
	
		
	}
	public void addRootTag(String sgName){
		for(int i=0; i<subGroupNames.size(); i++){
			if(subGroupNames.get(i).equals(sgName)){
				subGroupNames.set(i, "Root"+subGroupNames.get(i));
			}
		}
	}
	public void setPeaks(List<Point> ps){
		peaks.addAll(ps);
		for(Point p : peaks){
			regions.add(p.expand(win/2));
		}
	}
	public void setArffOutFileName(String arffOut){outArffFileName = arffOut;}
	public void setKmin(int kmin){Kmin = kmin;}
	public void setKmax(int kmax){Kmax = kmax;}
	public void setWin(int w){win = w;}
	
	//Getters
	public List<Region> getRandRegs(int n){
		List<Region> ret = new ArrayList<Region>();
		RandRegionsGenerator randReger = new RandRegionsGenerator(true, n);
		ret.addAll(randReger.execute());
		return ret;
	}
	
	public void execute() throws Exception{
		// First count K-mers and generate the .mat file
		KmerCounter counter = new KmerCounter();
		counter.printPeakSeqKmerRange(Kmin, Kmax);
		
		// Now read the mat file and generate the arff file
		generateArff();
		
		// Also, print the desgin file
		BufferedWriter bw = new BufferedWriter(new FileWriter("SeqUnwinder.design"));
		bw.write(design);
		bw.close();
	}
	
	public void generateArff() throws Exception{
		
		//
		CSVLoader loader = new CSVLoader();
		// Set options
		loader.setNominalAttributes("last");
		loader.setStringAttributes("");
		loader.setMissingValue("?");
		loader.setFieldSeparator("\t");
		StringBuilder sb = new StringBuilder();
		for(String s : subGroupNames){
			sb.append(s);sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		loader.setNominalLabelSpecs(new String[] {sb.toString()});
		loader.setFile(new File("tmpCounts.mat"));
		Instances data = loader.getDataSet();
		
		//Set subgroup index
		if(data.classIndex() == -1)
			data.setClassIndex(data.numAttributes()-1);
		
		//First, get weight index
		int wInd = data.numAttributes()-2;
		// Now set weights
		for(int i=0; i<data.numInstances(); i++){
			double weight = data.instance(i).value(wInd);
			data.instance(i).setWeight(weight);
		}
		// Now delete the weight attribute
		data.deleteAttributeAt(wInd);
		
		//Save the arff file
		ArffSaver saver = new ArffSaver();
		saver.setFile(new File(outArffFileName));
		saver.setInstances(data);
		saver.writeBatch();
		
		
	}
	
	
	
	public static void main(String[] args) throws Exception{
		
		// Reading genmome objects
		GenomeConfig gc = new GenomeConfig(args);
		MakeArff arffmaker = new MakeArff(gc);
		ArgParser ap = new ArgParser(args);
		
		// Get peak window size
		int win = Args.parseInteger(args, "win", 150);
		arffmaker.setWin(win);
		
		// Reading peaks files
		String eventsFile = ap.getKeyValue("peaks");
		
		FileReader fr = new FileReader(eventsFile);
		BufferedReader br = new BufferedReader(fr);
		String line;
		List<String> labels = new ArrayList<String>();
		while((line = br.readLine()) != null){
			if(!line.contains("#")){
				String[] pieces = line.split("\t");
				labels.add(pieces[1]);
			}
		}
		br.close();
		
		List<Point> peaks = new ArrayList<Point>();
		peaks.addAll(RegionFileUtilities.loadPeaksFromPeakFile(gc.getGenome(), eventsFile, -1));
		
		// Set peaks and regions
		arffmaker.setPeaks(peaks);
		
		// See if random regions/background regions are provide, if not generate them
		List<Region> ranregs = new ArrayList<Region>();
		
		if(ap.hasKey("randregs")){
			ranregs.addAll(RegionFileUtilities.loadRegionsFromPeakFile(gc.getGenome(), ap.getKeyValue("randregs"), win));
		}else{
			// First find how many rand regions are needed
			int numRand = Integer.MAX_VALUE; // This should be the size of the subgroup with minimum no of. instances
			Set<String> subgroupNames = new HashSet<String>();
			for(String s : labels){
				if(subgroupNames.add(s)){}
			}
			
			for(String s : subgroupNames){
				if(Collections.frequency(labels, s) < numRand)
					numRand = Collections.frequency(labels, s);
			}
			ranregs.addAll(arffmaker.getRandRegs(numRand));
		}
		
		// Now set labels, design file string, and ranRegs to peaks
		// Check if we have labels
		if(ap.hasKey("simple")){
			arffmaker.setLabelsSimple(labels, ranregs);
		}else{
			arffmaker.setLabels(labels, ranregs);
		}
		
		// Now get K-mer params
		int kmin = Args.parseInteger(args, "Kmin", 4);
		int kmax = Args.parseInteger(args, "Kmax", 6);
		arffmaker.setKmin(kmin);
		arffmaker.setKmax(kmax);
		
		// Name to write the arff file
		String arffOut = Args.parseString(args, "arffOut", "out.arff");
		if(!arffOut.endsWith(".arff"))
			arffOut = arffOut+".arff";
		
		arffmaker.setArffOutFileName(arffOut);
		
		// Now execute: Prints the SeqUnwinder design file and the arff file
		arffmaker.execute();
		
			
	}
	
	public class MapKeyComparator<X extends Number> implements Comparator<String>{
		HashMap<String,X> map;
		
		public MapKeyComparator(HashMap<String,X> m) {
			map = m;
		}
		
		public List<String> getKeyList(){
			List<String> ret = new ArrayList<String>();
			for(String s : map.keySet()){
				ret.add(s);
			}
			return ret;
			
		}

		@Override
		public int compare(String o1, String o2) {
			//return (map.get(o1).intValue()).compareTo(map.get(o2).intValue());
			return Integer.compare(map.get(o1).intValue(), map.get(o2).intValue());
		}
		
	}
	
	public class RandRegionsGenerator{
		
		private int numSamples = 1000;
		private int validSamples=0;
		
		private RepeatMaskedGenerator<Region> repMask;
		private double genomeSize=0;
		private long [] chromoSize;
		private String [] chromoNames;
		private int numChroms=0;
		private Random rand = new Random();
		private double repPropLimit=0.5;
		private boolean screenRepeats=false;
		
		//Settors
		public void setNum(int n){numSamples=n;}
		public void setScreenRepeats(boolean s){screenRepeats=s;}
		
		public RandRegionsGenerator(boolean screenReps, int num) {
			repMask = new RepeatMaskedGenerator<Region>(gcon.getGenome());
			setScreenRepeats(screenReps);
			setNum(num);
			
			
		}
		
		
		public List<Region> execute() {
			
			List<Region>regList = new ArrayList<Region>();
			// First see how big the genome is:
			chromoSize = new long[gcon.getGenome().getChromList().size()];
			chromoNames = new String[gcon.getGenome().getChromList().size()];
			Iterator<NamedRegion> chroms = new ChromRegionIterator(gcon.getGenome());
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				genomeSize += (double) currentChrom.getWidth();
				chromoSize[numChroms] = currentChrom.getWidth();
				chromoNames[numChroms] = currentChrom.getChrom();
				// System.out.println(chromoNames[numChroms]+"\t"+chromoSize[numChroms]);
				numChroms++;
			}// System.out.println(genomeSize);

			// Now, iteratively generate random positions and check if they are
			// valid
			while (validSamples < numSamples) {
				Region potential;
				long randPos = (long) (1 + (rand.nextDouble() * genomeSize));
				// find the chr
				boolean found = false;
				long total = 0;
				for (int c = 0; c < numChroms && !found; c++) {
					if (randPos < total + chromoSize[c]) {
						found = true;
						if (randPos + win < total + chromoSize[c]) {
							potential = new Region(gcon.getGenome(), chromoNames[c], (int) (randPos - total), (int) (randPos+ win - total - 1));
							boolean regionOK = true;

							// screen repeats
							if (screenRepeats) {
								// is this overlapping a repeat?
								double repLen = 0;
								Iterator<RepeatMaskedRegion> repItr = repMask.execute(potential);
								while (repItr.hasNext()) {
									RepeatMaskedRegion currRep = repItr.next();
									if (currRep.overlaps(potential)) {
										repLen += (double) currRep.getWidth();
									}
								}
								if (repLen / (double) potential.getWidth() > repPropLimit)
									regionOK = false;

								// Is the sequence free from N's?
								String potSeq = seqgen.execute(potential);
								if (potSeq.indexOf('N') >= 0) {
									regionOK = false;
								}
							}
							// Screen dupicates
							for (Region r : regList) {
								if (potential.overlaps(r))
									regionOK = false;
							}

							// Screen for any exclude regions provided
							if (regions.size() != 0) {
								for (Region ex : regions) {
									if (potential.overlaps(ex)) {
										regionOK = false;
									}
								}
							}

							if (regionOK) {
								validSamples++;
								regList.add(potential);
								System.out.println(potential.getChrom() + ":"
										+ potential.getStart() + "-"
										+ potential.getEnd());
							}
						}
					}
					total += chromoSize[c];
				}
			}
			return (regList);
		}
		
		
	}
	
	public class KmerCounter{
		
		protected String outfilename = "tmpCounts.mat";
		
		public KmerCounter() {}
		
		
		// Print the the kmer counts (for a given range of value k) in the sequences
		// for each peak
		public void printPeakSeqKmerRange(int kmin, int kmax) throws IOException {
			
			int numK = 0;
			for (int k = kmin; k <= kmax; k++) {
				numK = numK + (int) Math.pow(4, k);
			}
			int[] kmerCounts = new int[numK];
			
			FileWriter of = new FileWriter(outfilename);
			BufferedWriter bw = new BufferedWriter(of);
			
			StringBuilder sb = new StringBuilder();
			// Printing the header line
			for (int k = kmin; k <= kmax; k++) {
				int N = (int) Math.pow(4, k);
				for (int i = 0; i < N; i++)
					sb.append("\t"+RegionFileUtilities.int2seq(i, k));
			}
			sb.deleteCharAt(0);
			sb.append("\t"+"Weight"+"\t"+"Label"+"\n");
			//bw.write("Weight"+"\t"+"Label"+"\n");
			bw.write(sb.toString());

			//for (Region r : regs) {
			for(int r=0; r<regions.size(); r++){
				for (int i = 0; i < numK; i++)
					kmerCounts[i] = 0;

				String seq = seqgen.execute(regions.get(r)).toUpperCase();
				// Check if the sequence (seq) contains any N's if present ignore
				// them
				if (seq.contains("N"))
					continue;

				int ind = 0;
				for (int k = kmin; k <= kmax; k++) {
					for (int i = 0; i < (seq.length() - k + 1); i++) {
						String currK = seq.substring(i, i + k);
						String revCurrK = SequenceUtils.reverseComplement(currK);
						int currKInt = RegionFileUtilities.seq2int(currK);
						int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
						int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
						kmerCounts[ind + kmer]++;
					}
					ind = ind + (int) Math.pow(4, k);
				}
				sb = new StringBuilder();
				for (int i = 0; i < numK; i++)
					sb.append("\t" + kmerCounts[i]);
				sb.deleteCharAt(0);
				sb.append("\t"+subGroupWeights.get(subGroupNames.get(r))+"\t"+subGroupNames.get(r)+"\n");
				bw.write(sb.toString());
			}
			bw.close();
		}
		

	}

}
