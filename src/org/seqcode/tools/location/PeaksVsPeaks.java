package org.seqcode.tools.location;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;


public class PeaksVsPeaks {
	
	protected Genome gen;
	protected List<Point> peaksA;
	protected List<Region> regionsA;
	protected List<Point> peaksB;
	protected List<Region> regionsB;
	protected List<Point> peaksNeg;
	protected List<Region> regionsNeg;
	protected List<List<Point>> peaksBSet;
	protected List<List<Region>> regionsBSet;
	
	protected int overlapD=200; // min distance between the closest peaks to call them overlapping 
	protected int win=200; // for expanding a peak into a peak region (regions not currently used)
	
	protected int numRand=1000000;
	protected boolean monteCarloTesting=false;
	protected int numMCSets=1000;
	protected double VLARGEOVER=10000;
	
	public PeaksVsPeaks(Genome g) {
		gen = g;
	}
	
	
	
	// Settors
	public void setPeaksA(List<Point> a){peaksA = a;}
	public void setPeaksB(List<Point> b){peaksB = b;}
	public void setRegsA(List<Region> a){regionsA = a;}
	public void setRegsB(List<Region>  b){regionsB = b;}
	public void setOverlapD(int o){overlapD =o;}
	public void setMonteCarloTesting(int b){monteCarloTesting=true; numMCSets= b;}
	public void setWin(int w){win = w;}
	public void setNumRand(int r){numRand = r;}
	public void setPeaksBSet(List<List<Point>> pset){peaksBSet = pset;}
	public void setRegionsBSet(List<List<Region>> rset){regionsBSet = rset;}
	
	/**
	 * Given a list of peaks; prints the number of other peaks in that list
	 * that are within a distance of rad from each peak
	 * @param rad
	 */
	public void printClosePeaksInAList(int rad){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<Point, Integer> closebyPeaks = new HashMap<Point,Integer>();
		for(String chrom : peaksAbyChrom.keySet()){
			List<Point> currChrPeaks = peaksAbyChrom.get(chrom);
			Collections.sort(currChrPeaks);
			for(Point pa : currChrPeaks){
				int count = 0;
				for(Point pAprime : currChrPeaks){
					if(!pa.getLocationString().equals(pAprime.getLocationString())){
						if(pa.distance(pAprime) < rad)
							count++;
						else
							break;
					}
				}
				closebyPeaks.put(pa, count);
			}
		}

		StringBuilder sb = new StringBuilder();
		for(Point p : closebyPeaks.keySet()){
			sb.append(p.getLocationString());sb.append("\t");sb.append(closebyPeaks.get(p));sb.append("\n");
		}
		sb.deleteCharAt(sb.length()-1);
		System.out.println(sb.toString());
	}

	
	public void printUnique(int rad){
		
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		for(String chrom : peaksAbyChrom.keySet()){
			List<Point> outPeaks = new ArrayList<Point>();
			for(Point pa : peaksAbyChrom.get(chrom)){
				boolean add = true;
				for(Point addedA : outPeaks){
					if(addedA.expand(rad).contains(pa)){
						add = false;
						break;
					}
				}
				if(add)
					outPeaks.add(pa);
			}

			for(Point outP : outPeaks){
				System.out.println(outP.getLocationString());
			}

		}

	}
	
	public void printClosestOffsetPeak(){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<String, List<Point>> peaksBbyChrom = hashbychrom(peaksB);
		for(String chrom : peaksAbyChrom.keySet()){
			for(Point pa : peaksAbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				Point nearestPeak=null;
				for(Point pb : peaksBbyChrom.get(chrom)){
					if(pb.getChrom().equals(pa.getChrom())){
						if(pa.distance(pb) < Math.abs(mind)){
							mind = pa.offset(pb);
						nearestPeak = pb;
						}
					}
				}
				if(Math.abs(mind) <= overlapD){
					System.out.println(pa.getLocationString()+"\t"+nearestPeak.getLocationString()+"\t"+Integer.toString(mind));
				}else{
					System.out.println(pa.getLocationString()+"\tNA");
				}
			}

		}

	}
	
	
	public void printClosestPeaks(){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<String, List<Point>> peaksBbyChrom = hashbychrom(peaksB);
		for(String chrom : peaksAbyChrom.keySet()){
			if(!peaksBbyChrom.containsKey(chrom)){
				for(Point pa : peaksAbyChrom.get(chrom)){
					System.out.println(pa.getLocationString()+"\t-");
				}
				continue;
			}
				
			for(Point pa : peaksAbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				Point nearestPeak=null;
				for(Point pb : peaksBbyChrom.get(chrom)){
					if(pb.getChrom().equals(pa.getChrom())){
						if(pa.distance(pb) < mind){
							mind = pa.distance(pb);
							nearestPeak = pb;
						}
					}
				}
				if(mind <= overlapD){
					System.out.println(pa.getLocationString()+"\t"+nearestPeak.getLocationString()+"\t"+Integer.toString(mind));
				}else{
					System.out.println(pa.getLocationString()+"\t-");
				}
			
			}
		}
		
	}
	
	public void printOverlapStats(){
		System.out.println("PeaksATotal\tPeaksAOverlapB\tPeaksAOverlapBRate\tNegTotal\tNegOverlapB\tNegOverlapBRate\tPeaksAOverlapBOverRep");
		
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<String, List<Point>> peaksBbyChrom = hashbychrom(peaksB);
		HashMap<String, List<Point>> peaksNegbyChrom = hashbychrom(peaksNeg);
		//Counters
		double posPeaks=0, posPeaksRate;
		double negPeaks=0, negPeaksRate;
		double posTotal = (double)regionsA.size(), negTotal = (double)regionsNeg.size();
		
		//Positives
		for(String chrom : peaksAbyChrom.keySet()){
			for(Point pa : peaksAbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				if(peaksBbyChrom.containsKey(chrom)){
					for(Point pb : peaksBbyChrom.get(chrom)){
						if(pb.getChrom().equals(pa.getChrom())){
							if(pa.distance(pb) < mind){
								mind = pa.distance(pb);
							}
						}
					}
					if(mind <= overlapD)
						posPeaks++;
				}
			}
		}
		posPeaksRate = posPeaks/posTotal;
		
		//Negatives
		for(String chrom : peaksNegbyChrom.keySet()){
			for(Point pa : peaksNegbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				if(peaksBbyChrom.containsKey(chrom)){
					for(Point pb : peaksBbyChrom.get(chrom)){
						if(pb.getChrom().equals(pa.getChrom())){
							if(pa.distance(pb) < mind){
								mind = pa.distance(pb);
							}
						}
					}
					if(mind <= overlapD)
						negPeaks++;
				}
			}
		}
		negPeaksRate = negPeaks/negTotal;
		
		double peaksOverRep = negPeaksRate>0 ? posPeaksRate/negPeaksRate : VLARGEOVER;
		
		System.out.println(String.format("%.0f\t%.0f\t%.4f\t%.0f\t%.0f\t%.4f\t%.4f", posTotal,posPeaks,posPeaksRate,negTotal,negPeaks,negPeaksRate,peaksOverRep));
		
	}
	
	public void printOverlapStatsMC(){
		System.out.println("PeaksATotal\tPeaksAOverlapB\tPeaksAOverlapBRate\tNegMeanOverlapB\tNegMeanOverlapBRate\tNegStdDevOverlapBRate\tPeaksAOverlapBOverRep\tPeaksAOverlapBZScore");
		
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<String, List<Point>> peaksBbyChrom = hashbychrom(peaksB);
		HashMap<String, List<Point>> peaksNegbyChrom = hashbychrom(peaksNeg);
		//Counters
		double posPeaks=0, posPeaksRate;
		double negTotalPeaks=0, negPeaks=0, negPeaksRate, negStdDevPeaks=0, negMeanPeaksRate=0;
		double posTotal = (double)regionsA.size(), negTotal = (double)regionsNeg.size();
		
		//Positives
		for(String chrom : peaksAbyChrom.keySet()){
			for(Point pa : peaksAbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				if(peaksBbyChrom.containsKey(chrom)){
					for(Point pb : peaksBbyChrom.get(chrom)){
						if(pb.getChrom().equals(pa.getChrom())){
							if(pa.distance(pb) < mind){
								mind = pa.distance(pb);
							}
						}
					}
					if(mind <= overlapD)
						posPeaks++;
				}
			}
		}
		posPeaksRate = posPeaks/posTotal;
		
		//Negatives
		List<Double> negPeakRates = new ArrayList<Double>(); 
		int currStart;
		for(int n=0; n<numMCSets; n++){
			currStart=n*peaksA.size();
			negPeaks=0;
			List<Point> peaksNegSub = new ArrayList<Point>();
			for(int s=currStart; s<currStart+peaksA.size(); s++)
				peaksNegSub.add(peaksNeg.get(s));
			HashMap<String, List<Point>> peaksNegSubbyChrom = hashbychrom(peaksNegSub);
			
			for(String chrom : peaksNegSubbyChrom.keySet()){
				for(Point pa : peaksNegSubbyChrom.get(chrom)){
					int mind=Integer.MAX_VALUE;
					if(peaksBbyChrom.containsKey(chrom)){
						for(Point pb : peaksBbyChrom.get(chrom)){
							if(pb.getChrom().equals(pa.getChrom())){
								if(pa.distance(pb) < mind){
									mind = pa.distance(pb);
								}
							}
						}
						if(mind <= overlapD)
							negPeaks++;
					}
				}
			}
			negTotalPeaks+=negPeaks;
			negPeaksRate = negPeaks/posTotal;
			negPeakRates.add(negPeaksRate);
			negPeaksRate = negPeaks/negTotal;
		}
		
		//Negative stats
		double negMeanPeaks=negTotalPeaks/(double)(numMCSets);
		negMeanPeaksRate=negMeanPeaks/posTotal;
		double sqDiffSum=0;
		for(Double x : negPeakRates){
			sqDiffSum += (x-negMeanPeaksRate)*(x-negMeanPeaksRate);
		}negStdDevPeaks=Math.sqrt(sqDiffSum/(double)numMCSets);
		double peaksZScore = (posPeaks-negMeanPeaks)/negStdDevPeaks;
		double peaksOverRep = negMeanPeaksRate>0 ? posPeaksRate/negMeanPeaksRate : VLARGEOVER;
		
		System.out.println(String.format("%.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f", posTotal,posPeaks,posPeaksRate,negMeanPeaks,negMeanPeaksRate,negStdDevPeaks,peaksOverRep,peaksZScore));
		
	}
	
	public void printPeakSets(){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		List<HashMap<String,List<Point>>> peaksBSetbyChrom  = new ArrayList<HashMap<String,List<Point>>>();
		for(List<Point> set : peaksBSet){
			peaksBSetbyChrom.add(hashbychrom(set));
		}
		
		for(String chrom : peaksAbyChrom.keySet()){
			for(Point pa: peaksAbyChrom.get(chrom)){
				StringBuilder sb= new StringBuilder();
				sb.append(pa.getLocationString());sb.append("\t");
				Point[] nearestPeaks = new Point[peaksBSet.size()];
				int[] nearestDis = new int[peaksBSet.size()];
				Arrays.fill(nearestDis, Integer.MAX_VALUE);
				int setId=0;
				for(HashMap<String,List<Point>> bset: peaksBSetbyChrom){
					if(bset.containsKey(chrom)){
						for(Point pb : bset.get(chrom)){
							if(pb.getChrom().equals(pa.getChrom())){
								if(pa.distance(pb) < nearestDis[setId]){
									nearestDis[setId] = pa.distance(pb);
									nearestPeaks[setId] =pb;
								}
							}
						}
					}
					setId++;
				}
				for(int i=0; i<nearestDis.length;i++){
					if(nearestDis[i]<=overlapD){
						sb.append(nearestPeaks[i].getLocationString());sb.append("\t");
					}else{
						sb.append("-");sb.append("\t");
					}
				}
				sb.deleteCharAt(sb.length()-1);
				System.out.println(sb.toString());
			}
		}
	}

	public HashMap<String, List<Point>> hashbychrom(List<Point> pts){
		HashMap<String, List<Point>> byChr = new HashMap<String, List<Point>>();
		for(Point p : pts){
			if(!byChr.containsKey(p.getChrom()))
				byChr.put(p.getChrom(), new ArrayList<Point>());
			byChr.get(p.getChrom()).add(p);
		}
		return byChr;
	}
	
	//load negative
	public void loadNegative(String name){
		if(name==null || name.equals("random")){
			if(monteCarloTesting){
				numRand=peaksA.size()*numMCSets;
			}
			regionsNeg = randomRegionPick(regionsA, numRand, win);			
		}else{
			regionsNeg = RegionFileUtilities.loadRegionsFromPeakFile(gen, name, win);

			if(monteCarloTesting){
				numMCSets = regionsNeg.size()/regionsA.size();
				//System.out.println(numMCSets+" sets loaded for MonteCarlo testing");
			}
		}
		peaksNeg = new ArrayList<Point>();
		for(Region r : regionsNeg)
			peaksNeg.add(r.getMidpoint());
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
	
	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		if(!ap.hasKey("species") || !ap.hasKey("peaksA")|| !ap.hasKey("peaksB")) { 
            System.err.println("Usage:\n " +
                               "PeaksVsPeaks \n" +
                               " Required: \n" +
                               "  --species <species;genome> \n" +
                               "  --peaksA <file containing coordinates of peaks> or --peaksAset <file containing list of files> \n" +
                               "  --peaksB <file containing coordinates of peaks> or --peaksBset <file containing list of files> \n" +
                               "  --neg <random/filename> \n" +
                               "  --montecarlo <test this number of negative sets>\n" +
                               " More Info: \n"+
                               "  --overlapdist <d> \n"+
                               "  --numrand <number of random sequences to sample (ignored if using montecarlo option)> \n" +
                               " Options: \n" +
                               "  --peaksoverlap [print A peaks that overlap B peaks (default behavior)] \n" +
                               "  --printoffset [print the offset distance between overlapping peaks]\n" +
                               "  --overlapstats [freq/overrep statistics] \n" +
                               "  --uniq [print unique peaks in list A within --rad distance]\n"+
                               "  --rad <radius to use for some options>\n" +
                               "");
            return;
        }
		
		GenomeConfig gconfig = new GenomeConfig(args);
		
		int win = Args.parseInteger(args, "win", 150);
		
		PeaksVsPeaks analyzer = new PeaksVsPeaks(gconfig.getGenome());
		
		String peaksAfile = Args.parseString(args, "peaksA", null);
		if(peaksAfile == null){
			System.err.println("Provide peaks file A!!");
			return;
		}
		List<Point> peaksA = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		List<Region> regionsA = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		analyzer.setPeaksA(peaksA);
		analyzer.setRegsA(regionsA);
		
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksBset") && !ap.hasKey("uniq") && !ap.hasKey("closebyCount")){
			System.err.println("Provide peaks file to compare with the A list!!");
			return;
		}
		
		String peaksBfile = null;
		if(ap.hasKey("peaksB")){
			peaksBfile = Args.parseString(args, "peaksB", null);
			List<Point> peaksB = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksBfile, win);
			List<Region> regionsB = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksBfile, win);
			analyzer.setPeaksB(peaksB);
			analyzer.setRegsB(regionsB);
		}
		
		Collection<String> peaksBfileList = null;
		if(ap.hasKey("peaksBset")){
			peaksBfileList =  Args.parseStrings(args, "peaksBset");
			List<List<Point>> bPeaksSet = new ArrayList<List<Point>>();
			List<List<Region>> bRegionSet = new ArrayList<List<Region>>();
			for(String s : peaksBfileList){
				List<Point> currPoints = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), s);
				bPeaksSet.add(currPoints);
				List<Region> currRegions = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), s, win);
				bRegionSet.add(currRegions);
			}
			analyzer.setPeaksBSet(bPeaksSet);
			analyzer.setRegionsBSet(bRegionSet);
		}
		
		//Other options
		
		int overlapD = Args.parseInteger(args, "overlapdist", 200);
		analyzer.setOverlapD(overlapD);
		
		String neg = ap.hasKey("neg") ? ap.getKeyValue("neg"):null;
		if(ap.hasKey("overlapstats") && neg==null)
			neg="random";
		int numRand = Args.parseInteger(args, "numrand", 100000);
		analyzer.setNumRand(numRand);
		boolean runMC=false;
		if(ap.hasKey("montecarlo")){
			int mcTest = Args.parseInteger(args, "montecarlo", 1000);
			analyzer.setMonteCarloTesting(mcTest);
			runMC=true;
		}
		
		//Load Negatives if necessary
		if(neg!=null){
			analyzer.loadNegative(neg);
		}
		
		if(ap.hasKey("peaksB") && (ap.hasKey("peaksoverlap") || (!ap.hasKey("printoffset") && !ap.hasKey("overlapstats"))))
			analyzer.printClosestPeaks();
		else if(ap.hasKey("peaksB") && (ap.hasKey("overlapstats"))){
			if(runMC)
				analyzer.printOverlapStatsMC();
			else
				analyzer.printOverlapStats();
		}else if(ap.hasKey("peaksBset"))
			analyzer.printPeakSets();
		
		if(ap.hasKey("peaksB") && ap.hasKey("printoffset"))
			analyzer.printClosestOffsetPeak();
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksBset") && ap.hasKey("uniq")){
			int rad = Args.parseInteger(args, "rad", 50);
			analyzer.printUnique(rad);
		}
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksBset") && ap.hasKey("closebyCount")){
			int rad = Args.parseInteger(args, "rad", 50);
			analyzer.printClosePeaksInAList(rad);
		}
	}
	
	
	
	

}
