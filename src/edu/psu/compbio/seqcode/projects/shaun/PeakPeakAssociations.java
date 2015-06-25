package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import cern.jet.random.Binomial;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.projects.gps.features.EnrichedFeature;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
/**
 * Related to ReaksVsPeaks, for Peaks Vs Peaks statistics. However, we will keep this code static and limited 
 * in scope to allow reproduction of the RAR paper comparisons. 
 */
public class PeakPeakAssociations {

	private Species org=null;
	private Genome gen =null;
	private ArrayList<Region> posSet;
	private ArrayList<Point> posPeaks;
	private ArrayList<Region> negSet;
	private HashMap<String,ArrayList<Region>> testSet = new HashMap<String,ArrayList<Region>>();
	private HashMap<String,ArrayList<Point>> testPeaks = new HashMap<String,ArrayList<Point>>();
	private int numRand=1000000;
	private int window=200;
	private boolean monteCarloTesting=false;
	private int numMCSets=500;
	private final double VLARGEOVER=3000000000.0;
	private int histBinSize=10;
	private double [] histogram=null;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")|| !ap.hasKey("peaks")|| (!ap.hasKey("testpeaks")&&!ap.hasKey("testdomains"))) { 
            System.err.println("Usage:\n " +
                               "PeakPeakAssociations \n" +
                               " Required: \n" +
                               "  --species <organism name> " +
                               "  --genome <genome version> "+
                               "  --peaks <file containing coordinates of primary peaks> \n" +
                               "  --neg <random/filename> \n" +
                               "  --montecarlo [flag for MonteCarlo testing]\n" +
                               "  --nummc <number of MC sets>" +
                               "  --printneg [flag to print negative regions to file]" +
                               "  --testpeaks <peaks to compare to> \n" +
                               "  --testdomains <domains to compare to> \n" +
                               " More Info: \n"+
                               "  --win <window of sequence around positive/negative points> \n"+
                               "  --numrand <number of random sequences to sample> \n" +
                               " Options: \n" +
                               "  --overlapstats [freq/overrep statistics for each peak] \n" +
                               "  --bitpattern [print present/absent bit patterns] \n" +
                               "  --printhisto [print peak-peak distance histogram] \n" +
                               "");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
    	String posFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks"):null;
        String neg = ap.hasKey("neg") ? ap.getKeyValue("neg"):null;
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        int numSamp = 1000000;
        if(ap.hasKey("numrand")){
        	numSamp = new Integer(ap.getKeyValue("numrand")).intValue();
        }
        Collection<String> testPeaks = ap.hasKey("testpeaks") ? Args.parseStrings(args, "testpeaks") : null;
        Collection<String> testDomains = ap.hasKey("testdomains") ? Args.parseStrings(args, "testdomains") : null;
        
        //options
        boolean peaksOverlap = ap.hasKey("peaksoverlap");
        boolean overlapStats = ap.hasKey("overlapstats");
        boolean bitPattern = ap.hasKey("bitpattern");
        boolean mcTest = ap.hasKey("montecarlo");
        boolean printNeg = ap.hasKey("printneg");
        boolean printDistHisto = ap.hasKey("printhisto");
        
        
        try {
			Species currorg = Species.getSpecies(species);
			Genome currgen = new Genome(currorg, genome);
        
			//initialize
			PeakPeakAssociations analyzer = new PeakPeakAssociations(currorg, currgen);
			
			//load options
			analyzer.setMonteCarloTesting(mcTest);
			if(ap.hasKey("nummc")){
	        	analyzer.setNumMC(new Integer(ap.getKeyValue("nummc")).intValue());
	        }
			analyzer.setNumTest(numSamp);
			analyzer.setWin(win);
			
			//load positive & negative sets
			analyzer.loadPositive(posFile);
			analyzer.loadNegative(neg);
			analyzer.loadTest(testPeaks, testDomains);

			//Options
			if(printNeg)
				analyzer.printNegativeRegions();
			if(printDistHisto)
				analyzer.printDistHistogram();
			if(overlapStats){
				if(analyzer.isMonteCarloTesting())
					analyzer.printOverlapStatsMC();
				else
					analyzer.printOverlapStats();
			}if(bitPattern)
				analyzer.printBitPattern();
			
			//analyzer.printMotifInfo();			
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public PeakPeakAssociations(Species o, Genome g){
		org = o;
		gen=g;
	}
	///////////////////////////////////////////////////////////////////////
	//Options first
	///////////////////////////////////////////////////////////////////////


	//Histogram of peak-peak distances
	public void printDistHistogram(){
		int numBins = window/histBinSize;
		histogram = new double [numBins+1];
		for(int b=0; b<=numBins; b++){ histogram[b]=0;}
		
		for(String t : testSet.keySet()){
			ArrayList<Point> currTest = testPeaks.get(t); 
			
			for(int s=0; s<posPeaks.size(); s++){
				Point p = posPeaks.get(s);
				for(Point pt : currTest)
					if(p.getChrom().equals(pt.getChrom()) && p.distance(pt)<=window){
						int dist = p.distance(pt);
						int bin = dist/histBinSize;
						histogram[bin]++;
					}
			}
		}
		for(int b=0; b<=numBins; b++){
			int bin = b*histBinSize;
			System.out.println(bin+"\t"+histogram[b]);
		}
	}
	//Print a bit patterns for each peak
	public void printBitPattern(){
		boolean [][] contains = new boolean[posSet.size()][testSet.keySet().size()];
		for(int i=0; i<posSet.size(); i++)
			for(int j=0; j<testSet.keySet().size(); j++)
				contains[i][j]=false;
		
		int y=0;
		for(String t : testSet.keySet()){
			ArrayList<Region> currTest = testSet.get(t); 
			
			for(int s=0; s<posSet.size(); s++){
				Region r = posSet.get(s);
				for(int x=0; x<currTest.size(); x++){
					Region rt = currTest.get(x);
					if(r.overlaps(rt))
						contains[s][y]=true;
				}
			}
			y++;
		}
		
		System.out.print("Expt");
		for(String t : testSet.keySet()){System.out.print("\t"+t);}
		System.out.print("\n");
		for(int i=0; i<posSet.size(); i++){
			System.out.print(posSet.get(i));
			int z=0;
			for(String t : testSet.keySet()){
				if(contains[i][z])
					System.out.print("\t1");
				else
					System.out.print("\t0");
				z++;
			}System.out.print("\n");
		}
	}
	//Print some occurrence and over-representation info for each set
	public void printOverlapStats(){
		ArrayList<Result> results = new ArrayList<Result>();
		System.out.println("Expt\tPosTotal\tPosPeaks\tPosPeaksRate\tNegTotal\tNegPeaks\tNegPeaksRate\tPeakOverRep\tPVal");
		double posTotal = (double)posSet.size(), negTotal = (double)negSet.size();
		for(String t : testSet.keySet()){
			ArrayList<Region> currTestSet = testSet.get(t);
			
			//Counters
			double posPeaks=0, posPeaksRate;
			double negPeaks=0, negPeaksRate;
			
			//Positive set
			for(int s=0; s<posSet.size(); s++){
				Region r = posSet.get(s);
				boolean goodPeak =false;
				for(Region q : currTestSet){
					if(r.overlaps(q))
						goodPeak=true;
				}if(goodPeak)
					posPeaks++;
			}posPeaksRate = posPeaks/posTotal;
			
			//Negative set
			for(int s=0; s<negSet.size(); s++){
				Region r = negSet.get(s);
				boolean goodPeak =false;
				for(Region q : currTestSet){
					if(r.overlaps(q))
						goodPeak=true;
				}if(goodPeak)
					negPeaks++;
			}negPeaksRate = negPeaks/negTotal;
			
			double peaksOverRep = negPeaksRate>0 ? posPeaksRate/negPeaksRate : VLARGEOVER;
			
			double posnegScale = negTotal/posTotal;
			double PVal;
			if(negPeaks==0 && posPeaks==0)
				PVal=1;
			else
				PVal =  binomialPValue(negPeaks/posnegScale,posPeaks+(negPeaks/posnegScale));
						
			String desc = new String(t+"\t"+posTotal+"\t"+posPeaks+"\t"+posPeaksRate+"\t"+negTotal+"\t"+"\t"+negPeaks+"\t"+negPeaksRate+"\t"+peaksOverRep+"\t");
			results.add(new Result(desc, PVal));
			//System.out.println(t+"\t"+posTotal+"\t"+posPeaks+"\t"+posPeaksRate+"\t"+negTotal+"\t"+"\t"+negPeaks+"\t"+negPeaksRate+"\t"+peaksOverRep+"\t"+PVal);
		}
		Collections.sort(results);
		results = benjaminiHochbergCorrection(results);
		Collections.sort(results);
		
		for(Result r : results)
			System.out.println(r.desc+String.format("%e", r.pval));
	}
	//Print some occurrence and over-representation info for each set
	public void printOverlapStatsMC(){
		System.out.println("Expt\tExptTotal\tPosTotal\tPosPeaks\tPosPeaksRate\tNegMeanPeaks\tNegMeanPeaksRate\tNegStdDevPeaksRate\tPeakOverRep\tP-Value");
		double posTotal = (double)posSet.size(), negTotal = (double)negSet.size(), negOneTotal =(double)negSet.size()/(double)numMCSets; 
		for(String t : testSet.keySet()){
			ArrayList<Region> currTestSet = testSet.get(t);
			
			//Counters
			double peaksZScore=0;
			double posPeaks=0, posPeaksRate;
			double negPeaks=0, negPeaksRate=0, negTotalPeaks=0, negStdDevPeaks=0,negVarPeaks, negMeanPeaksRate=0;
			
			//Positive set
			for(int s=0; s<posSet.size(); s++){
				Region r = posSet.get(s);
				boolean goodPeak =false;
				for(Region q : currTestSet){
					if(r.overlaps(q)){
						goodPeak=true;
						break;
					}
				}if(goodPeak)
					posPeaks++;
			}posPeaksRate = posPeaks/posTotal;
			
			//Negative set
			ArrayList<Double> negPeakRates = new ArrayList<Double>(); 
			int currStart;
			for(int n=0; n<numMCSets; n++){
				currStart=n*posSet.size();
				negPeaks=0;
				for(int s=currStart; s<currStart+negOneTotal; s++){
					Region r = negSet.get(s);
					boolean goodPeak =false;
					for(Region q : currTestSet){
						if(r.overlaps(q)){
							goodPeak=true;
							break;
						}
					}if(goodPeak)
						negPeaks++;
				}negTotalPeaks+=negPeaks;
				negPeaksRate = negPeaks/negOneTotal;
				negPeakRates.add(negPeaksRate);
			}
			
			//Negative stats
			double negMeanPeaks=negTotalPeaks/(double)(numMCSets);
			negMeanPeaksRate=negMeanPeaks/negOneTotal;
			double sqDiffSum=0;
			for(Double x : negPeakRates){
				sqDiffSum += (x-negMeanPeaksRate)*(x-negMeanPeaksRate);
			}negVarPeaks = sqDiffSum/(double)numMCSets;
			negStdDevPeaks=Math.sqrt(sqDiffSum/(double)numMCSets);
			peaksZScore = (posPeaks-negMeanPeaks)/negStdDevPeaks;
			double peaksOverRep = negMeanPeaksRate>0 ? posPeaksRate/negMeanPeaksRate : VLARGEOVER;
			
			Normal norm = new Normal(negMeanPeaks, negVarPeaks, new DRand());
			double PVal = 1-norm.cdf(posPeaks);
			
			System.out.println(String.format("%s\t%d\t%.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%e", t,currTestSet.size(),posTotal,posPeaks,posPeaksRate,negMeanPeaks,negMeanPeaksRate,negStdDevPeaks,peaksOverRep,PVal));
			//System.out.println(t+"\t"+currTestSet.size()+"\t"+posTotal+"\t"+posPeaks+"\t"+posPeaksRate+"\t"+negMeanPeaks+"\t"+negMeanPeaksRate+"\t"+negStdDevPeaks+"\t"+peaksOverRep+"\t"+peaksZScore);
		}
	}
	
	///////////////////////////////////////////////////////////////////////
	public void setMonteCarloTesting(boolean mc){monteCarloTesting=mc;}
	public boolean isMonteCarloTesting(){return monteCarloTesting;}
	public void setNumTest(int n){numRand=n;}
	public void setNumMC(int mc){numMCSets=mc;}
	public void setWin(int w){window=w;}

	//load positive
	public void loadPositive(String fname){
		posSet = loadRegionsFromPeakFile(fname, window);
		posPeaks = loadPeaksFromPeakFile(fname, window);
	}
	//load negative
	public void loadNegative(String name){
		if(name==null || name.equals("random")){
			if(isMonteCarloTesting()){
				numRand=posSet.size()*numMCSets;
			}
			negSet = randomRegionPick(posSet, numRand, window);			
		}else{
			if(isMonteCarloTesting()){
				numRand=posSet.size()*numMCSets;
				ArrayList<Region> tmpSet = loadRegionsFromPeakFile(name, window);
				negSet = randomRegionSample(tmpSet, numRand);
			}else{
				negSet = loadRegionsFromPeakFile(name, window);
			}
		}
	}
	//load test
	public void loadTest(Collection<String> peakF, Collection<String> domainF){
		if(peakF!=null){
			for(String f : peakF){
				String n = f.replaceAll(".peaks", "");
				String[] x = n.split("/");
				testSet.put(x[x.length-1], loadRegionsFromPeakFile(f, window));
				testPeaks.put(x[x.length-1], loadPeaksFromPeakFile(f, window));
			}
		}
		if(domainF !=null){
			for(String f : domainF){
				String n = f.replaceAll(".domains", "");
				String[] x = n.split("/");
				testSet.put(x[x.length-1], loadRegionsFromPeakFile(f, -1));
				testPeaks.put(x[x.length-1], loadPeaksFromPeakFile(f, window));
			}
		}
	}
	
	//Load a set of regions from a peak file
	public ArrayList<Region> loadRegionsFromPeakFile(String filename, int win){
		ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(!words[0].equals("Region")){
		            if(words.length>=3 && win!=-1){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
	                	int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
	                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
	                	regs.add(r);
	                }else if(words.length>=1){
		            	RegionParser parser = new RegionParser(gen);
		            	Region r = parser.execute(words[0]);
		            	if(r!=null){regs.add(r);}
		            }
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(regs);
	}
	//Load a set of peaks from a peak file
	public ArrayList<Point> loadPeaksFromPeakFile(String filename, int win){
		ArrayList<Point> points = new ArrayList<Point>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(!words[0].equals("Region")){
		            if(words.length>=3 && win!=-1){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	points.add(p);
	                }else if(words.length>=1){
		            	RegionParser parser = new RegionParser(gen);
		            	Region r = parser.execute(words[0]);
		            	if(r!=null){points.add(new Point(r.getGenome(), r.getChrom(), (r.getStart()+r.getEnd())/2));}
		            }
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(points);
	}

	//Randomly pick a set of Regions
	private ArrayList<Region> randomRegionPick(ArrayList<Region> blackList, int numSamples, int sampleSize){
		ArrayList<Region> regs = new ArrayList<Region>();
		Random rand = new Random();
		int validSamples=0;
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[gen.getChromList().size()];
		String [] chromoNames = new String[gen.getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		while(validSamples<numSamples){
			Region potential;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			boolean found=false;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					if(randPos+sampleSize<total+chromoSize[c]){
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total));
						
						//is this region in the blacklist? 
						boolean valid=true;
						if(blackList!=null){
							for(Region r : blackList){
								if(potential.overlaps(r)){valid=false;}
							}
						}
						if(valid){
							validSamples++;
							regs.add(potential);
						}
					}
				}total+=chromoSize[c];
			}
		}
		return(regs);
	}
	//Randomly sample from a defined set of regions
	private ArrayList<Region> randomRegionSample(ArrayList<Region> fullSet, int num){
		ArrayList<Region> samples = new ArrayList<Region>();
		Random rand = new Random();
		for(int i=0; i<num; i++){
			int r = rand.nextInt(fullSet.size());
			samples.add(fullSet.get(r));
		}
		return(samples);
	}
	
	public void printNegativeRegions(){
		FileWriter fout = null;
		try {
			fout = new FileWriter("neg.regions");

			for(Region n : negSet){
				fout.write(n.getLocationString()+"\n");
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		finally {
			if (fout != null) {
				try {
					fout.close();
				}
				catch (IOException ioex) {
					ioex.printStackTrace();
				}
			}
		}
	}
	
	//Binomial CDF assuming scaled control. Uses COLT binomial test
	// k=scaled control, n=scaled control+signal
	protected double binomialPValue(double k, double n){
		double pval=1;
		Binomial b = new Binomial((int)Math.ceil(n), 0.5, new DRand());
		pval = b.cdf((int) Math.ceil(k));
		return(pval);		
	}
	//Multiple hypothesis testing correction -- assumes results ordered according to p-value
	protected ArrayList<Result> benjaminiHochbergCorrection(ArrayList<Result> res){
		double total = res.size();
		ArrayList<Result> out = new ArrayList<Result>();
		double rank =1;
		for(Result r : res){
			r.pval = r.pval*(total/rank);
			if(r.pval>1)
				r.pval=1;
			out.add(r);
			rank++;
		}return(res);
	}
	
	private class Result implements Comparable<Result>{
		public String desc;
		public double pval;
		
		public Result(String d, double p){
			desc=d;
			pval=p;
		}
		//Rank according to increasing p-value 
		public int compareTo(Result r) {
			if(pval<r.pval){return(-1);}
			else if(pval>r.pval){return(1);}
			else{return(0);}
		}
	}
}
