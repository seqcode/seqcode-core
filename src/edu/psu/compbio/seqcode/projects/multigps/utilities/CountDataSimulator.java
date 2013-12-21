package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.stats.NegativeBinomialDistrib;

/**
 * CountDataSimulator: generate a two-condition CountsDataset file given an empirically observed dataset and some parameters
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CountDataSimulator {

	//Parameters
	private int numReplicates=1;
	private final int numConditions=2;
	private int numDataPoints = 10000;
	private Double alpha = 0.15; //Parameter for biological variance. Poisson distribution (technical variance only) is alpha=0
	private Double upRegFrac=0.1; //Fraction of genes that are up-regulated B vs A
	private Double downRegFrac=0.1; //Fraction of genes that are down-regulated B vs A
	private Double diffExpLevel =2.0; //Base fold-change of differentially regulated genes
	private Double condATotalReads=10000000.0;
	private Double condBTotalReads=20000000.0;
	
	//Variables
	private ArrayList<Double> empirical = new ArrayList<Double>();
	
	//Output
	double[][] counts;
	double[][] extraCounts;
	List<SimCounts> simResults = new ArrayList<SimCounts>();
	
	//Constructor
	public CountDataSimulator(){}
	
	//Modifiers
	public void setDataPoints(int nd){numDataPoints=nd;}
	public void setReplicates(int r){numReplicates=r;}
	//public void setConditions(int c){numConditions=c;}
	public void setAlpha(Double a){alpha=a;}
	public void setReadsA(Double r){condATotalReads = r;}
	public void setReadsB(Double r){condBTotalReads = r;}
	public void setUpRegFrac(Double f){upRegFrac=f;}
	public void setDownRegFrac(Double f){downRegFrac=f;}
	public void setDiffExpLevel(Double d){diffExpLevel = d;}
	
	/**
	 * Load an empirical dataset from a file
	 * @param filename
	 */
	public void loadEmpiricalFromFile(String filename){
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line= reader.readLine(); //Skip first line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            if(!line.startsWith("#")){
	            	String[] words = line.split("\\s+");
	            	double e = new Double(words[1]);
	            	empirical.add(e);
	            }
	        }reader.close();
	        Collections.sort(empirical);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Simulate counts using a negative binomial model
	 * @param outFilename
	 */
	public List<SimCounts> simulate(){
		counts = new double[numDataPoints][numConditions*numReplicates];
		extraCounts = new double[numDataPoints][numConditions*numReplicates];
		double[] absolutes = new double[numDataPoints];
		simResults = new ArrayList<SimCounts>();
		boolean [] diffs = new boolean[numDataPoints]; 
		Random rand = new Random();
		NegativeBinomialDistrib nb = new NegativeBinomialDistrib(1,0.5);
		double [] condAMols = new double[numDataPoints];
		double [] condBMols = new double[numDataPoints];
		double [] condMolTotals = new double[numConditions];
		double [] sampleReadTotals = new double[numConditions*numReplicates];
		double [] sampleReadsPerMols = new double[numConditions*numReplicates];
		double absTotal = 0;
		int x=0;
		for(int c=0; c<numConditions; c++){
			condMolTotals[c]=0;
			for(int r=0; r<numReplicates; r++){
				sampleReadTotals[x]=0;
				sampleReadsPerMols[x]=0;
				x++;
			}
		}
		
		//Assume that empirical holds the absolute molecule counts
		//Generate new molecule counts for two conditions
		for(int d=0; d<numDataPoints; d++){
			//Choose a binding event strength
			double dice = rand.nextDouble();
			int strengthIndex = (int)(dice*(double)empirical.size());
			double empMol = empirical.get(strengthIndex);
			absolutes[d] = empMol;

			double currAMol = empMol;
			double currBMol = empMol;
			dice = rand.nextDouble();
			//Is this gene differentially expressed in this condition? (w.r.t condition A)
			//Implemented a little strangely to avoid skew
			diffs[d] = true;
			if(dice<(downRegFrac/2))
				currAMol=currAMol*diffExpLevel;
			else if(dice<downRegFrac)
				currBMol=currBMol/diffExpLevel;
			else if(dice<(downRegFrac+upRegFrac/2))
				currAMol=currAMol/diffExpLevel;
			else if(dice<(downRegFrac+upRegFrac))
				currBMol=currBMol*diffExpLevel;
			else
				diffs[d] = false;
			
			condAMols[d]=currAMol;
			condBMols[d]=currBMol;
			condMolTotals[0]+=currAMol;
			condMolTotals[1]+=currBMol;
			absTotal+=empMol;
		}

		//For each condition
		for(int c=0; c<numConditions; c++){	
			//Sample counts at the appropriate concentrations
			for(int d=0; d<numDataPoints; d++){
				double conc = c==0 ? condAMols[d]/condMolTotals[0] : condBMols[d]/condMolTotals[1];
				double mean = c==0 ? (conc*condATotalReads) : (conc*condBTotalReads);
				double var = mean+(mean*mean*alpha);
				nb.setMeanVar(mean, var);
				int sample = c*numReplicates;
				for(int r=0; r<numReplicates; r++){
					counts[d][sample]=(double)nb.nextInt();
					sampleReadTotals[sample]+=counts[d][sample];
					sample++;
				}
			}
		}
		
		//Extra counts around the absolute conc
		for(int d=0; d<numDataPoints; d++){
			double conc = absolutes[d]/absTotal;
			for(int c=0; c<numConditions; c++){
				double mean = c==0 ? (conc*condATotalReads) : (conc*condBTotalReads);
				double var = mean+(mean*mean*alpha);
				nb.setMeanVar(mean, var);
				int sample= c*numReplicates;
				for(int r=0; r<numReplicates; r++){
					extraCounts[d][sample]=(double)nb.nextInt();
					sample++;
				}
			}
		}
		
		//Make simulated results data structure
		for(int d=0; d<numDataPoints; d++){
			SimCounts sim = new SimCounts();
			sim.absolute = absolutes[d];
			sim.counts=counts[d];
			sim.backup = extraCounts[d];
			sim.isDiff=diffs[d];
			simResults.add(sim);
		}
		
		//Reads per Mol ratios
		System.out.println("CountDataSimulator: Reads per Mol ratios\n");
		x=0;
		for(int c=0; c<numConditions; c++){
			for(int r=0; r<numReplicates; r++){
				System.out.print("\t"+c+":"+r);
				sampleReadsPerMols[x]=sampleReadTotals[x]/condMolTotals[c];
				x++;
			}
		}System.out.println("");
		x=0;
		for(int c1=0; c1<numConditions; c1++){
			for(int r1=0; r1<numReplicates; r1++){
				int y=0;
				System.out.print(c1+":"+r1);
				for(int c2=0; c2<numConditions; c2++){
					for(int r2=0; r2<numReplicates; r2++){
						System.out.print(String.format("\t%.2f",(sampleReadsPerMols[x]/sampleReadsPerMols[y])));
						y++;
					}
				}System.out.println("");
				x++;
			}
		}System.out.println("");
		
		return(simResults);
	}
	
	/**
	 * Write the counts table to a file
	 * @param outFilename
	 */
	public void printOutput(String outBase, boolean printLabels){
		if(counts!=null){
			try {
				String outFilename = outBase+".counts";
				FileWriter fw = new FileWriter(outFilename);
				fw.write("Point");
				if(printLabels)
					fw.write("\tLabel");
				for(int c=0; c<numConditions; c++){
					for(int r=0; r<numReplicates; r++)
						fw.write("\t"+c+":"+r);
				}fw.write("\n");
				for(int d=0; d<numDataPoints; d++){
					fw.write(d+"");
					if(printLabels){
						int label = simResults.get(d).isDiff? 1:0;
						fw.write("\t"+label);
					}
					for(int s=0; s<(numConditions*numReplicates); s++)
						fw.write("\t"+counts[d][s]);
					fw.write("\n");
				}
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}else{
			System.err.println("Counts matrix not defined.");
		}
	}
	
	/**
	 * Main file for simulating counts data from the command-line
	 * @param args
	 */
	public static void main(String[] args) {
		String empFile, outFile="out";
		int r, numdata;
		double rA, rB, a, up, down, diff;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h") || !ap.hasKey("emp")){
			System.err.println("CountDataSimulator:\n" +
					"\t--emp <empirical data file>\n" +
					"\t--numdata <number of data points>\n" +
					"\t--r <num replicates per condition>\n" +
					"\t--a <over-dispersion param>\n" +
					"\t--readsA <avg reads in cond A>\n" +
					"\t--readsB <avg reads in cond B>\n" +
					"\t--up <up-regulated fraction>\n" +
					"\t--down <down-regulated fraction>\n" +
					"\t--diff <basis of differential expression>\n" +
					"\t--out <output file>\n" +
					"");
		}else{
			CountDataSimulator sim = new CountDataSimulator();
			
			if(ap.hasKey("emp")){
				empFile = ap.getKeyValue("emp");
				sim.loadEmpiricalFromFile(empFile);
			}if(ap.hasKey("numdata")){
				numdata = new Integer(ap.getKeyValue("numdata"));
				sim.setDataPoints(numdata);
			}if(ap.hasKey("r")){
				r = new Integer(ap.getKeyValue("r"));
				sim.setReplicates(r);
			}if(ap.hasKey("a")){
				a = new Double(ap.getKeyValue("a"));
				sim.setAlpha(a);
			}if(ap.hasKey("readsA")){
				rA = new Double(ap.getKeyValue("readsA"));
				sim.setReadsA(rA);
			}if(ap.hasKey("readsB")){
				rB = new Double(ap.getKeyValue("readsB"));
				sim.setReadsB(rB);
			}if(ap.hasKey("up")){
				up = new Double(ap.getKeyValue("up"));
				sim.setUpRegFrac(up);
			}if(ap.hasKey("down")){
				down = new Double(ap.getKeyValue("down"));
				sim.setDownRegFrac(down);
			}if(ap.hasKey("diff")){
				diff = new Double(ap.getKeyValue("diff"));
				sim.setDiffExpLevel(diff);
			}if(ap.hasKey("out")){
				outFile = ap.getKeyValue("out");
			}
			
			sim.simulate();
			sim.printOutput(outFile, false);
		}
	}
	
	public class SimCounts{
		double absolute; //The absolute value of binding enrichment (mean before differential)
		double [] counts; //The simulated counts
		double [] backup; //Extra counts used by the multi-condition read simulator
		boolean isDiff=false;
	}
}


