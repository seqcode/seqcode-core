package edu.psu.compbio.seqcode.projects.multigps.utilities; 

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.ReadHit;
import edu.psu.compbio.seqcode.projects.multigps.utilities.CountDataSimulator.SimCounts;

/**
 * Simulates multi-condition reads using BindingModels. <br> 
 * 
 * In any given simulated replicate experiment, reads are simulated according to a BindingModel.
 * Read counts for each binding event in each condition are determined from a CountDataSimulator,
 * which allows for differential binding between two conditions, and negative-binomial sampled
 * read counts in each replicate. 
 * Noise reads are distributed uniformly across the genome. 
 * 
 * @author Shaun Mahony
 *
 */
public class MultiConditionReadSimulator {

	
	private BindingModel model;
	private int numConditions = 2;
	private int numReplicates = 1;
	private String outPath;
	private FileWriter[][] writers;
	
	private Genome fakeGen;
	private long[] chromLens;
	private HashMap<String, Long> chromOffsets = new HashMap<String, Long>();
	private int rLen=32;
	private long genomeLength=-1;
	private double noiseProbabilities[][];
	private int numSigReads[][], numTotalReads[][];
	private List<SimCounts> simCounts;
	
	private int numEvents=0;
	private int eventSpacing = 10000;
	private double jointEventRate = 0.0;
	private int jointEventSpacing = 200;
	private List<Pair<Point, SimCounts>> events = new ArrayList<Pair<Point, SimCounts>>();
	private HashMap<Point, Boolean> eventIsJoint = new HashMap<Point, Boolean>();
	
	public MultiConditionReadSimulator(BindingModel m, Genome g, List<SimCounts> counts, int numCond, int numRep, double noiseProb, double jointRate, int jointSpacing, String outPath){
		model=m;
		numConditions = numCond;
		numReplicates = numRep;
		jointEventRate = jointRate;
		jointEventSpacing = jointSpacing;
		this.outPath = outPath;
		fakeGen = g;
		genomeLength = (long)fakeGen.getGenomeLength();
		chromLens = new long[fakeGen.getChromList().size()];
		
		int c=0; long offset=0;
		for(String chr : fakeGen.getChromList()){
			chromLens[c]=fakeGen.getChromLength(chr); 
			chromOffsets.put(chr, offset);
			offset+=chromLens[c];
			c++;
		}
		
		noiseProbabilities = new double[numConditions][numReplicates];
		numSigReads = new int[numConditions][numReplicates];
		numTotalReads = new int[numConditions][numReplicates];
		for(int co=0; co<numConditions; co++)
			for(int r=0; r<numReplicates; r++){
				noiseProbabilities[co][r]=noiseProb;
				numSigReads[co][r]=0;
				numTotalReads[co][r]=0;
			}
		
		try {
			writers = new FileWriter[numConditions][numReplicates];
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++){
					FileWriter fout = new FileWriter(outPath+"_reads_C"+co+"_R"+r+".bed");
					writers[co][r]=fout;
				}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if(noiseProb<1.0){
			simCounts = counts;
			setBindingPositions();
			numEvents = events.size();
			
			//Count signal reads and infer total reads
			for(Pair<Point,SimCounts> event : events){
				int index=0;
				for(int co=0; co<numConditions; co++)
					for(int r=0; r<numReplicates; r++){
						if(!eventIsJoint.get(event.car())){
							numSigReads[co][r]+=event.cdr().counts[index];
						}else{
							numSigReads[co][r]+=event.cdr().backup[index];
						}
						numTotalReads[co][r] = (int)((double)numSigReads[co][r]/(1-noiseProbabilities[co][r]));
						index++;
					}
			}
		}
	}
	
	/**
	 * Make binding positions corresponding to all simulated counts.
	 * By default, shared events are in chr1, and diff events are in chrX
	 */
	private void setBindingPositions(){
		Random jointDice = new Random();
		long sharedOffset=chromOffsets.get("1"), diffOffset=chromOffsets.get("X");
		
		
		for(SimCounts s : simCounts){
			long curroff =0;
			if(s.isDiff){
				curroff =diffOffset;
				diffOffset+=eventSpacing;
			}else{
				curroff =sharedOffset;
				sharedOffset+=eventSpacing;
			}
			
			//Translate from offsets to chromosome name and start
			int c=0; long offset=0;
			String chr = fakeGen.getChromList().get(0);
			while(curroff>(offset+chromLens[c]) && c<fakeGen.getChromList().size()-1){
				c++;
				chr = fakeGen.getChromList().get(c);
				offset = chromOffsets.get(chr);
			}
			long start = curroff-offset;
			Point p = new Point(fakeGen, chr, (int)start);
			events.add(new Pair<Point, SimCounts>(p,s));
			eventIsJoint.put(p, false);
			
			//Simulate a joint event?
			double jointrand = jointDice.nextDouble();
			if(jointrand<jointEventRate){
				Point jp = new Point(fakeGen, chr, (int)start+jointEventSpacing);
				events.add(new Pair<Point, SimCounts>(jp,s));
				eventIsJoint.put(jp, true);
			}
		}
	}
	
	/**
	 * Simulate a set of noise reads for each replicate in each condition.
	 * Numbers of reads are predefined;
	 */
	private void simulateNoiseReads(){
		try {
			Random noiseGenerator = new Random();
			Random strandGenerator = new Random();
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++){
					List<ReadHit> reads = new ArrayList<ReadHit>();
					int noiseReads = (int)((double)numTotalReads[co][r]*noiseProbabilities[co][r]);
					for(int i=0; i<noiseReads; i++){
						ReadHit rh=null;
						double noiserand = noiseGenerator.nextDouble();
						double strandrand = strandGenerator.nextDouble();
	
						long pos = (long)(noiserand*(genomeLength));
						
						//Translate from pos to chromosome name and start
						int c=0; long offset=0;
						String chr = fakeGen.getChromList().get(0);
						while(pos>(offset+chromLens[c]) && c<fakeGen.getChromList().size()-1){
							c++;
							chr = fakeGen.getChromList().get(c);
							offset = chromOffsets.get(chr);
						}
						long start = pos-offset;
						
						//Add the ReadHit
						if (strandrand<0.5)
							rh = new ReadHit(chr, (int)start, (int)start+rLen-1, '+');
						else
							rh = new ReadHit(chr, Math.max(1, (int)start-rLen+1), (int)start, '-');
						reads.add(rh);
					}
					for(ReadHit rh : reads){
						writers[co][r].write(rh.getChrom()+"\t"+rh.getStart()+"\t"+rh.getEnd()+"\tU\t0\t"+rh.getStrand()+"\n");
					}
				}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Simulate a set of binding event reads for each replicate in each condition. 
	 * Number of reads, number of events, and strengths of events are all pre-determined. 
	 */
	private void simulateBindingReads(){
		//Initialize the probability landscape
		int eventWidth=1000; int evoff = eventWidth/2;
		double[] forProbLand=new double[eventWidth]; double[] revProbLand=new double[eventWidth];
		double[] forProbCumul=new double[eventWidth]; double[] revProbCumul=new double[eventWidth];
		for(int i=0; i<eventWidth; i++){
			forProbLand[i]=0; revProbLand[i]=0;
			forProbCumul[i]=0; revProbCumul[i]=0;
		}
		int modelRange = Math.max(Math.abs(model.getMin()), Math.abs(model.getMax()));
		int winStart = Math.max(0, evoff-modelRange);
		int winStop = Math.min(evoff+modelRange, eventWidth-1);
		for(int i=winStart; i<winStop; i++){
			int forDist  = i-evoff;
			int revDist  = evoff-i;
			forProbLand[i]+=model.probability(forDist);
			revProbLand[i]+=model.probability(revDist);
		}
		//Set the cumulative scores
		double fTotal=0, rTotal=0;
		for(int i=0; i<eventWidth; i++){
			fTotal+=forProbLand[i];
			rTotal+=revProbLand[i];
			forProbCumul[i]=fTotal;
			revProbCumul[i]=rTotal;
		}
		//Normalize
		for(int i=0; i<eventWidth; i++){
			forProbLand[i]=forProbLand[i]/fTotal;
			revProbLand[i]=revProbLand[i]/rTotal;
			forProbCumul[i]=forProbCumul[i]/fTotal;
			revProbCumul[i]=revProbCumul[i]/rTotal;
		}
		
		
		//Generate the reads
		try {
			Random sigGenerator = new Random();
			Random strandGenerator = new Random();
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++){
					List<ReadHit> reads = new ArrayList<ReadHit>();
					
					// Generate reads for each event
					for(Pair<Point, SimCounts> ps : events){
						Point coord = ps.car();
						String chr = coord.getChrom();
						int chrLen = fakeGen.getChromLength(chr);
						SimCounts sc = ps.cdr();
						int sample = co*numReplicates+r;
						double readCount = eventIsJoint.get(coord) ? sc.backup[sample] : sc.counts[sample];
						
						ReadHit rh = null;
						for(int x=0; x<readCount; x++){
							boolean forwardStrand = strandGenerator.nextDouble()<0.5 ? true:false;
							double rand = sigGenerator.nextDouble();
							int fivePrimeEnd=coord.getLocation()-evoff;
							//Find the probability interval
							if(forwardStrand){
								for(int j=0; j<eventWidth; j++){
									if(forProbCumul[j] > rand){
										fivePrimeEnd+=j;
										break;
									}
								}
								//Make the ReadHit
								if(fivePrimeEnd<chrLen)
									rh = new ReadHit(chr, fivePrimeEnd, fivePrimeEnd+rLen-1, '+');			
							}else{
								for(int j=eventWidth-1; j>=0; j--){
									if(revProbCumul[j] < rand){
										fivePrimeEnd+=j+1;
										break;
									}
								}
								//Make the ReadHit
								rh = new ReadHit(chr, Math.max(1, fivePrimeEnd-rLen+1), fivePrimeEnd, '-');
							}
							reads.add(rh);
						}
					}
					//Write the reads to the file
					for(ReadHit rh : reads){
						writers[co][r].write(rh.getChrom()+"\t"+rh.getStart()+"\t"+rh.getEnd()+"\tU\t0\t"+rh.getStrand()+"\n");
					}
				}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
		
	//Accessors
	public void setTotalReads(int totReadsA, int totReadsB){
		for(int co=0; co<numConditions; co++)
			for(int r=0; r<numReplicates; r++){
				if(co==0)
					numTotalReads[co][r] = totReadsA;
				else
					numTotalReads[co][r] = totReadsB;
			}
	}
	public void setJointEventRate(double jointRate){
		this.jointEventRate = jointRate;
	}

	// clean up
	public void close(){
		try {
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++){
					writers[co][r].close();
				}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Write an output file of event positions & expected read counts per replicate
	 */
	public void printEvents(){
		try {
			String filename = outPath+".events";
			FileWriter fout = new FileWriter(filename);
			
			fout.write("Coord\tDiff");
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++)
					fout.write("\tC"+co+"_R"+r);
			fout.write("\n");
			
			// Generate reads for each event
			for(Pair<Point, SimCounts> ps : events){
				Point coord = ps.car();
				SimCounts sc = ps.cdr();
				int diff = sc.isDiff && !eventIsJoint.get(coord) ? 1 : 0;
				fout.write(coord.getLocationString()+"\t"+diff);
				for(int x=0; x<sc.counts.length; x++){
					if(!eventIsJoint.get(coord))
						fout.write("\t"+sc.counts[x]);
					else
						fout.write("\t"+sc.backup[x]);
				}fout.write("\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**

	 */
	public static void main(String[] args) {
		String empFile, outFile="out";
		int r=2, numdata, jointSpacing=200;
		double rA=1000000, rB=100000, a, up, down, diff, jointRate=0.0;
		String bmfile;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h") || !ap.hasKey("emp")){
			System.err.println("MultiConditionReadSimulator:\n" +
					"\t--geninfo <genome info file>\n" +
					"\t--emp <empirical data file>\n" +
					"\t--numdata <number of events to simulate>\n" +
					"\t--r <num replicates per condition>\n" +
					"\t--a <over-dispersion param>\n" +
					"\t--readsA <avg total reads in cond A>\n" +
					"\t--readsB <avg total reads in cond B>\n" +
					"\t--up <up-regulated fraction>\n" +
					"\t--down <down-regulated fraction>\n" +
					"\t--diff <basis of differential expression>\n" +
					"\t--model <binding model file>\n" +
					"\t--noise <noise probability per replicate>\n" +
					"\t--jointrate <proportion of peaks that are joint binding events>\n" +
					"\t--jointspacing <spacing between joint events>\n" +
					"\t--out <output file>\n" +
					"");
		}else{
			Genome gen = null;
			CountDataSimulator cdsim = new CountDataSimulator();
			
			//////////////////////////////////////////////////
			// Read in parameters 
			if(Args.parseArgs(args).contains("geninfo")){
				//Make fake genome... chr lengths provided
				String fName = Args.parseString(args, "geninfo", null);
				gen = new Genome("Genome", new File(fName), true);
			}
			if(ap.hasKey("emp")){
				empFile = ap.getKeyValue("emp");
				cdsim.loadEmpiricalFromFile(empFile);
			}if(ap.hasKey("numdata")){
				numdata = new Integer(ap.getKeyValue("numdata"));
				cdsim.setDataPoints(numdata);
			}if(ap.hasKey("r")){
				r = new Integer(ap.getKeyValue("r"));
				cdsim.setReplicates(r);
			}if(ap.hasKey("a")){
				a = new Double(ap.getKeyValue("a"));
				cdsim.setAlpha(a);
			}if(ap.hasKey("readsA")){
				rA = new Double(ap.getKeyValue("readsA"));
			}if(ap.hasKey("readsB")){
				rB = new Double(ap.getKeyValue("readsB"));
				cdsim.setReadsB(rB);
			}if(ap.hasKey("up")){
				up = new Double(ap.getKeyValue("up"));
				cdsim.setUpRegFrac(up);
			}if(ap.hasKey("down")){
				down = new Double(ap.getKeyValue("down"));
				cdsim.setDownRegFrac(down);
			}if(ap.hasKey("diff")){
				diff = new Double(ap.getKeyValue("diff"));
				cdsim.setDiffExpLevel(diff);
			}if(ap.hasKey("jointrate")){
				jointRate = new Double(ap.getKeyValue("jointrate"));
			}if(ap.hasKey("jointspacing")){
				jointSpacing = new Integer(ap.getKeyValue("jointspacing"));
			}if(ap.hasKey("out")){
				outFile = ap.getKeyValue("out");
			}
			double noiseProb   = Args.parseDouble(args, "noise", 0.9);
			double [][] noiseProbs = new double[2][r];
			for(int x=0; x<2; x++)
				for(int y=0; y<noiseProbs.length; y++)
					noiseProbs[x][y]=noiseProb;
			cdsim.setReadsA(rA*(1-noiseProb));
			cdsim.setReadsB(rB*(1-noiseProb));
			
			bmfile  = Args.parseString(args, "model", null);
			File mFile = new File(bmfile);
			if(!mFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
			
			
			//////////////////////////////////////////////////
			// Simulate counts 
			List<SimCounts> counts = null;
			if(noiseProb<1.0){
				counts = cdsim.simulate();
				cdsim.printOutput(outFile, true);
			}
			
			//////////////////////////////////////////////////
			// Simulate reads according to counts and binding model
			BindingModel bm = new BindingModel(mFile);
	        //Initialize the MultiConditionReadSimulator
	        MultiConditionReadSimulator sim = new MultiConditionReadSimulator(bm, gen, counts, 2, r, noiseProb, jointRate, jointSpacing, outFile);
	        if(noiseProb==1.0)
	        	sim.setTotalReads((int)rA, (int)rB);

	        if(noiseProb<1){
	        	sim.printEvents();
	        	sim.simulateBindingReads();
	        }
			sim.simulateNoiseReads();
			sim.close();
		}
	}
}
