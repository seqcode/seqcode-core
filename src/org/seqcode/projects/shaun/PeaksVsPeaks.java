package org.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gse.gsebricks.verbs.location.PointParser;
import org.seqcode.gse.gsebricks.verbs.location.RegionParser;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;
import org.seqcode.gse.utils.io.RegionFileUtilities;


public class PeaksVsPeaks {

	private Species org=null;
	private Genome gen =null;
	private List<Region> posSet;
	private List<Point> posPeaks;
	private List<Region> negSet;
	private HashMap<String,List<Region>> testSet = new HashMap<String,List<Region>>();
	private HashMap<String,List<Point>> testPeaks = new HashMap<String,List<Point>>();
	private List<String> posLines;
	private int numRand=1000000;
	private int window=200;
	private boolean monteCarloTesting=false;
	private int numMCSets=1000;
	private final double VLARGEOVER=3000000000.0;
	private int histBinSize=10;
	private double [] histogram=null;
	private boolean absolute=false;
	private boolean strandedRegions=false;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("peaks")|| (!ap.hasKey("testpeaks")&&!ap.hasKey("testdomains"))) { 
            System.err.println("Usage:\n " +
                               "PeaksVsPeaks \n" +
                               " Required: \n" +
                               "  --species <species;genome> " +
                               "  --peaks <file containing coordinates of primary peaks> \n" +
                               "  --neg <random/filename> \n" +
                               "  --montecarlo [flag for MonteCarlo testing]\n" +
                               "  --printneg [flag to print negative regions to file]" +
                               "  --testpeaks <peaks to compare to> \n" +
                               "  --testdomains <domains to compare to> \n" +
                               " More Info: \n"+
                               "  --win <window of sequence around positive/negative points> \n"+
                               "  --numrand <number of random sequences to sample> \n" +
                               " Options: \n" +
                               "  --peaksoverlap [peaks containing ANY other peak] \n" +
                               "  --overlapstats [freq/overrep statistics for each peak] \n" +
                               "  --bitpattern [print present/absent bit patterns] \n" +
                               "  --printhisto [print peak-peak distance histogram] \n" +
                               "  --printclosesthisto [print peak-peak closest distance histogram] \n" +
                               "  --printclosestpeakpairs [print the peak and closest test peak if in window] \n" +
                               "  --binsize <bin size for the histogram> \n" +
                               "  --absolute [flag for absolute distances]\n" +
                               "  --sr [flag for stranded region format]\n" +
                               "");
            return;
        }
        
        try {
	        Pair<Species, Genome> pair = Args.parseGenome(args);
	    	Species currorg = pair.car();
	    	Genome currgen = pair.cdr();
	    	String posFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks"):null;
	        String neg = ap.hasKey("neg") ? ap.getKeyValue("neg"):null;
	        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
	        int binSize = ap.hasKey("binsize") ? new Integer(ap.getKeyValue("binsize")).intValue():10;
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
	        boolean printClosestDistHisto = ap.hasKey("printclosesthisto");
	        boolean printClosestPeakPairs = ap.hasKey("printclosestpeakpairs");
	        boolean abs = ap.hasKey("absolute");
	        boolean sr = ap.hasKey("sr");
	        
			//initialize
			PeaksVsPeaks analyzer = new PeaksVsPeaks(currorg, currgen);
			
			//load options
			analyzer.setMonteCarloTesting(mcTest);
			analyzer.setNumTest(numSamp);
			analyzer.setWin(win);
			analyzer.setAbsolute(abs);
			analyzer.setStrandedRegions(sr);
			analyzer.setBinSize(binSize);
			
			//load positive & negative sets
			analyzer.loadPositive(posFile);
			analyzer.loadNegative(neg);
			analyzer.loadTest(testPeaks, testDomains);

			//Options
			if(printNeg)
				analyzer.printNegativeRegions();
			if(peaksOverlap)
				analyzer.printPeaksOverlap();
			if(printDistHisto)
				analyzer.printDistHistogram();
			if(printClosestDistHisto)
				analyzer.printClosestDistHistogram();
			if(printClosestPeakPairs)
				analyzer.printClosestPeakPairs();
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
	
	public PeaksVsPeaks(Species o, Genome g){
		org = o;
		gen=g;
	}
	///////////////////////////////////////////////////////////////////////
	//Options first
	///////////////////////////////////////////////////////////////////////

	//Simple printing of peak lines that contain ANY of the peaks
	public void printPeaksOverlap(){
		boolean [] contains = new boolean[posSet.size()];
		for(int i=0; i<posSet.size(); i++){contains[i]=false;}
		for(String t : testSet.keySet()){
			List<Region> currTest = testSet.get(t); 
			
			for(int s=0; s<posSet.size(); s++){
				Region r = posSet.get(s);
				for(Region rt : currTest)
					if(r.overlaps(rt))
						contains[s]=true;
			}
		}
		for(int i=0; i<posSet.size(); i++){
			if(contains[i])
				System.out.println(posLines.get(i));
		}
	}
	//Histogram of peak-peak distances
	public void printDistHistogram(){
		int numBins = (window/histBinSize)*2;
		int zeroBin = (numBins/2);
		histogram = new double [numBins+1];
		for(int b=0; b<=numBins; b++){ histogram[b]=0;}
		
		for(String t : testSet.keySet()){
			List<Point> currTest = testPeaks.get(t); 
			HashMap<String, List<Point>> currByChr = hashPeaksByChr(currTest);
			
			for(int s=0; s<posPeaks.size(); s++){
				Point p = posPeaks.get(s);
				if(currByChr.containsKey(p.getChrom())){
					for(Point pt : currByChr.get(p.getChrom()))
						if(p.distance(pt)<=window){
							int dist = p.distance(pt);
							int bin=zeroBin;
							if(!absolute){
								if(p.getLocation()<=pt.getLocation()) 
									bin = zeroBin+(dist/histBinSize);
								else
									bin = zeroBin-(dist/histBinSize);
							}else{
								bin = zeroBin+(dist/histBinSize);
							}
							histogram[bin]++;
						}
				}
			}
		}
		for(int b=0; b<zeroBin; b++){
			int bin = (zeroBin-b)*histBinSize;
			System.out.println(bin+"\t"+histogram[b]);
		}
		for(int b=zeroBin; b<=numBins; b++){
			int bin = (b-zeroBin)*histBinSize;
			System.out.println(bin+"\t"+histogram[b]);
		}
	}
	//Histogram of closest peak-peak distances
	public void printClosestDistHistogram(){
		int numBins = (window/histBinSize)*2;
		int zeroBin = (numBins/2);
		histogram = new double [numBins+1];
		for(int b=0; b<=numBins; b++){ histogram[b]=0;}
		
		for(String t : testSet.keySet()){
			List<Point> currTest = testPeaks.get(t); 
			HashMap<String, List<Point>> currByChr = hashPeaksByChr(currTest);
			
			//Find the nearest peak 
			for(int s=0; s<posPeaks.size(); s++){
				int closestDist = Integer.MAX_VALUE;
				Point closestPoint=null;
				Point p = posPeaks.get(s);
				if(currByChr.containsKey(p.getChrom())){
					for(Point pt : currByChr.get(p.getChrom())){
						if(Math.abs(p.distance(pt))<=closestDist){
							closestDist = Math.abs(p.distance(pt));
							closestPoint = pt;
						}
					}
					if(closestDist<=window && closestPoint !=null){
						int bin=zeroBin;
						if(!absolute){
							if(p.getLocation()<=closestPoint.getLocation()) 
								bin = zeroBin+(closestDist/histBinSize);
							else
								bin = zeroBin-(closestDist/histBinSize);
						}else{
							bin = zeroBin+(closestDist/histBinSize);
						}
						histogram[bin]++;
					}
				}
			}
		}
		for(int b=0; b<zeroBin; b++){
			int bin = (zeroBin-b)*histBinSize;
			System.out.println(bin+"\t"+histogram[b]);
		}
		for(int b=zeroBin; b<=numBins; b++){
			int bin = (b-zeroBin)*histBinSize;
			System.out.println(bin+"\t"+histogram[b]);
		}
	}
	//Print closest peak pairs with distance
	public void printClosestPeakPairs(){
		
		for(String t : testSet.keySet()){
			List<Point> currTest = testPeaks.get(t); 
			HashMap<String, List<Point>> currByChr = hashPeaksByChr(currTest);
			
			//Find the nearest peak 
			for(int s=0; s<posPeaks.size(); s++){
				int closestDist = Integer.MAX_VALUE;
				Point closestPoint=null;
				Point p = posPeaks.get(s);
				if(currByChr.containsKey(p.getChrom())){
					for(Point pt : currByChr.get(p.getChrom())){
						if(Math.abs(p.distance(pt))<=closestDist){
							closestDist = Math.abs(p.distance(pt));
							closestPoint = pt;
						}
					}
					if(closestDist<=window && closestPoint !=null){
						System.out.println(p.getLocationString()+"\t"+closestPoint.getLocationString()+"\t"+p.offset(closestPoint));
					}
				}
			}
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
			List<Region> currTest = testSet.get(t); 
			
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
	//Print some occurrence and over-representation info for each motif
	public void printOverlapStats(){
		System.out.println("Expt\tPosTotal\tPosPeaks\tPosPeaksRate\tNegTotal\tNegPeaks\tNegPeaksRate\tPeakOverRep");
		double posTotal = (double)posSet.size(), negTotal = (double)negSet.size();
		for(String t : testSet.keySet()){
			List<Region> currTestSet = testSet.get(t);
			
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
			
			System.out.println(t+"\t"+posTotal+"\t"+posPeaks+"\t"+posPeaksRate+"\t"+negTotal+"\t"+"\t"+negPeaks+"\t"+negPeaksRate+"\t"+peaksOverRep);
		}
	}
	//Print some occurrence and over-representation info for each motif
	public void printOverlapStatsMC(){
		System.out.println("Expt\tExptTotal\tPosTotal\tPosPeaks\tPosPeaksRate\tNegMeanPeaks\tNegMeanPeaksRate\tNegStdDevPeaksRate\tPeakOverRep\tPeakZScore");
		double posTotal = (double)posSet.size(), negTotal = (double)negSet.size();
		for(String t : testSet.keySet()){
			List<Region> currTestSet = testSet.get(t);
			
			//Counters
			double peaksZScore=0;
			double posPeaks=0, posPeaksRate;
			double negPeaks=0, negPeaksRate=0, negTotalPeaks=0, negStdDevPeaks=0, negMeanPeaksRate=0;
			
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
				for(int s=currStart; s<currStart+posSet.size(); s++){
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
				negPeaksRate = negPeaks/posTotal;
				negPeakRates.add(negPeaksRate);
			}
			
			//Negative stats
			double negMeanPeaks=negTotalPeaks/(double)(numMCSets);
			negMeanPeaksRate=negMeanPeaks/posTotal;
			double sqDiffSum=0;
			for(Double x : negPeakRates){
				sqDiffSum += (x-negMeanPeaksRate)*(x-negMeanPeaksRate);
			}negStdDevPeaks=Math.sqrt(sqDiffSum/(double)numMCSets);
			peaksZScore = (posPeaks-negMeanPeaks)/negStdDevPeaks;
			double peaksOverRep = negMeanPeaksRate>0 ? posPeaksRate/negMeanPeaksRate : VLARGEOVER;
			
			System.out.println(String.format("%s\t%d\t%.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f", t,currTestSet.size(),posTotal,posPeaks,posPeaksRate,negMeanPeaks,negMeanPeaksRate,negStdDevPeaks,peaksOverRep,peaksZScore));
			//System.out.println(t+"\t"+currTestSet.size()+"\t"+posTotal+"\t"+posPeaks+"\t"+posPeaksRate+"\t"+negMeanPeaks+"\t"+negMeanPeaksRate+"\t"+negStdDevPeaks+"\t"+peaksOverRep+"\t"+peaksZScore);
		}
	}
	
	///////////////////////////////////////////////////////////////////////
	public void setMonteCarloTesting(boolean mc){monteCarloTesting=mc;}
	public boolean isMonteCarloTesting(){return monteCarloTesting;}
	public void setNumTest(int n){numRand=n;}
	public void setWin(int w){window=w;}
	public void setBinSize(int b){histBinSize=b;}
	public void setAbsolute(boolean a){absolute=a;}
	public void setStrandedRegions(boolean s){strandedRegions=s;}

	//load positive
	public void loadPositive(String fname){
		if(strandedRegions){
			posSet= new ArrayList<Region>();
			posPeaks= new ArrayList<Point>();
			for(StrandedRegion sr : RegionFileUtilities.loadStrandedRegionsFromMotifFile(gen, fname, window))
				posSet.add(new Region(sr.getGenome(), sr.getChrom(), sr.getStart(), sr.getEnd()));
			for(StrandedPoint sp : RegionFileUtilities.loadStrandedPointsFromMotifFile(gen, fname, window))
				posPeaks.add(new Point(sp.getGenome(), sp.getChrom(), sp.getLocation()));
		}else{
			posSet = RegionFileUtilities.loadRegionsFromPeakFile(gen, fname, window);
			posPeaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, fname, window);
		}
		posLines = RegionFileUtilities.loadLinesFromFile(fname);
	}
	//load negative
	public void loadNegative(String name){
		if(name==null || name.equals("random")){
			if(isMonteCarloTesting()){
				numRand=posSet.size()*numMCSets;
				//System.out.println(numMCSets+" sets randomly generated for MonteCarlo testing");
			}
			negSet = randomRegionPick(posSet, numRand, window);			
		}else{
			if(strandedRegions){
				negSet = new ArrayList<Region>();
				for(StrandedRegion sr : RegionFileUtilities.loadStrandedRegionsFromMotifFile(gen, name, window))
					negSet.add(new Region(sr.getGenome(), sr.getChrom(), sr.getStart(), sr.getEnd()));
			}else{
				negSet = RegionFileUtilities.loadRegionsFromPeakFile(gen, name, window);
			}
			if(isMonteCarloTesting()){
				numMCSets = negSet.size()/posSet.size();
				System.out.println(numMCSets+" sets loaded for MonteCarlo testing");
			}
		}
	}
	//load test
	public void loadTest(Collection<String> peakF, Collection<String> domainF){
		if(peakF!=null){
			for(String f : peakF){
				String n = f.replaceAll(".peaks", "");
				String[] x = n.split("/");
				if(strandedRegions){
					List<Region> tmpSet = new ArrayList<Region>();
					List<Point> tmpPeaks = new ArrayList<Point>();
					for(StrandedRegion sr : RegionFileUtilities.loadStrandedRegionsFromMotifFile(gen, f, window))
						tmpSet.add(new Region(sr.getGenome(), sr.getChrom(), sr.getStart(), sr.getEnd()));
					for(StrandedPoint sp : RegionFileUtilities.loadStrandedPointsFromMotifFile(gen, f, window))
						tmpPeaks.add(new Point(sp.getGenome(), sp.getChrom(), sp.getLocation()));
					testSet.put(x[x.length-1], tmpSet);
					testPeaks.put(x[x.length-1], tmpPeaks);
				}else{
					testSet.put(x[x.length-1], RegionFileUtilities.loadRegionsFromPeakFile(gen, f, window));
					testPeaks.put(x[x.length-1], RegionFileUtilities.loadPeaksFromPeakFile(gen, f, window));
				}
			}
		}
		if(domainF !=null){
			for(String f : domainF){
				String n = f.replaceAll(".domains", "");
				String[] x = n.split("/");
				testSet.put(x[x.length-1], RegionFileUtilities.loadRegionsFromPeakFile(gen, f, -1));
				testPeaks.put(x[x.length-1], RegionFileUtilities.loadPeaksFromPeakFile(gen, f, window));
			}
		}
	}
	

	//Randomly pick a set of Regions
	private List<Region> randomRegionPick(List<Region> blackList, int numSamples, int sampleSize){
		List<Region> regs = new ArrayList<Region>();
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
	
	private HashMap<String, List<Point>> hashPeaksByChr(List<Point> peaks){
		HashMap<String, List<Point>> byChr = new HashMap<String, List<Point>>();
		for(Point p : peaks){
			if(!byChr.containsKey(p.getChrom()))
				byChr.put(p.getChrom(), new ArrayList<Point>());
			byChr.get(p.getChrom()).add(p);
		}
		return byChr;
	}
}
