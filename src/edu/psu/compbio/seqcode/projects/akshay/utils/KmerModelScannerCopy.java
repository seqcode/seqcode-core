package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.BayesmentsSandbox;
import edu.psu.compbio.seqcode.projects.akshay.clusterkmerprofile.ClusterProfiles;

public class KmerModelScannerCopy {
	
	protected GenomeConfig gcon;
	protected double[] kmerweights;
	protected double[][] kmerpairwights;
	protected boolean isPairModel = false;
	protected List<Point> posPeaks = new ArrayList<Point>();
	protected List<Region> posRegions = new ArrayList<Region>();
	protected List<Point> negPeaks = new ArrayList<Point>();
	protected List<Region> negRegions = new ArrayList<Region>();
	protected int minK; // length of the smallest K-mer
	protected int maxK; // length of the largest K-mer
	protected int minM; // Minimum length to consider for motif finding
	protected int maxM; //Maximum length to consider for motif finding
	//protected int n; // Number of motifs to look at each 
	protected String outbase;
	protected File outdir;
	
	protected List<Region> posMountains = new ArrayList<Region>();
	protected ArrayList<int[]> posProfiles = new ArrayList<int[]>();;
	protected HashMap<Integer,String> posMountainsToIndex =new HashMap<Integer,String>();
	protected HashMap<Integer,Double> posMountainsScores =new HashMap<Integer,Double>();
	
	protected List<Region> negMountains = new ArrayList<Region>();
	protected ArrayList<int[]> negProfiles =new ArrayList<int[]>();
	protected HashMap<Integer,Double> negMountainsScores = new HashMap<Integer,Double>();
	
	protected HashMap<Integer,String> negMountainsToIndex = new HashMap<Integer,String>();
	// All clustering parameters
	protected int its_CLUS=100;
	protected int numClus_CLUS=3;
	
	// Two options at the momment 
	// fixedwin :- A fixed window around all peaks 
	// sliding:- Sliding apprach with window size varying from mimM to maxM
	protected String SCAN_TYPE ="fixedwin";
	
	
	
	/** 
	 * Constructor
	 * @param gc
	 */
	public KmerModelScannerCopy(GenomeConfig gc) {
		gcon=gc;
	}
	
	
	
	// Settors
	public void setKmerWeights(double[] w){kmerweights=w;}
	public void setPosPeaks(List<Point>ps){posPeaks=ps;}
	public void setPosRegions(List<Region>rs){posRegions=rs;}
	public void setNegPeaks(List<Point>ps){negPeaks=ps;}
	public void setNegRegions(List<Region>rs){negRegions=rs;}
	public void setMinK(int K){minK=K;}
	public void setMaxK(int K){maxK=K;}
	public void setminM(int m){minM=m;}
	public void setmaxM(int m){maxM=m;}
	public void setModelType(boolean isPair){isPair = isPairModel;}
	public void setKmerPairWeights(double[][] pw){kmerpairwights = pw;}
	public void setOutbase(String o){outbase=o;}
	public void setOutDir(File odir){outdir=odir;}
	public void setClusterItrs(int numItrs){its_CLUS=numItrs;}
	public void setNumClusters(int numClus){numClus_CLUS = numClus;}
	public void setScanType(String scanType){SCAN_TYPE=scanType;}
	
	//Gettors
	public int getKmerBaseInd(String kmer){
		int len = kmer.length();
		int baseInd = 0;
		for(int k=minK; k<kmer.length(); k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	
	
	//Fillers
	private void fillMountains(boolean useCache, String genpath, double oddsThresh) throws IOException{
		List<Pair<Region,Double>> posMounts= new ArrayList<Pair<Region,Double>>();
		posMounts= findMountains(useCache,genpath,oddsThresh,posRegions);
		int index=0;
		HashMap<String,Integer> addedMountains = new HashMap<String,Integer>();
		for(Pair<Region,Double> pr : posMounts){
			if(!addedMountains.containsKey(pr.car().getLocationString())){
				posMountains.add(pr.car());
				posMountainsToIndex.put(index,pr.car().getLocationString() );
				posMountainsScores.put(index++, pr.cdr());
				addedMountains.put(pr.car().getLocationString(), 1);
			}
		}
		posProfiles = getProfilesAtPeaks(posMountains, useCache,genpath);
		
		index=0;
		addedMountains = new HashMap<String,Integer>();
		
		if(negPeaks.size() !=0 && negRegions.size() !=0){
			List<Pair<Region,Double>> negMounts = new ArrayList<Pair<Region,Double>>();
			negMounts= findMountains(useCache,genpath,-1*oddsThresh,negRegions);
		
			for(Pair<Region,Double> pr : negMounts){
				if(!addedMountains.containsKey(pr.car().getLocationString())){
					negMountains.add(pr.car());
					negMountainsToIndex.put(index,pr.car().getLocationString());
					negMountainsScores.put(index++, pr.cdr());
					addedMountains.put(pr.car().getLocationString(), 1);
				}
			}
			negProfiles =getProfilesAtPeaks(negMountains, useCache,genpath);
		}
		
	}
	
	
	
	/**
	 * Prints the mountain composition
	 * @param useCache
	 * @param genpath
	 * @param oddsThresh
	 */
	public void clusterKmerProfilesAtMountains(boolean useCache, String genpath, double oddsThresh) throws IOException{
		fillMountains(useCache,genpath,oddsThresh);
		//ClusterProfiles clusterManager = new ClusterProfiles(its_CLUS,numClus_CLUS,posProfiles,posMountainsToIndex,k,posMountainsScores,outbase,outdir);
		//clusterManager.execute("pos");
		//if(negPeaks.size() !=0 && negRegions.size() !=0){
		//	clusterManager = new ClusterProfiles(its_CLUS,numClus_CLUS,negProfiles,negMountainsToIndex,k,negMountainsScores,outbase,outdir);
		//	clusterManager.execute("neg");
		//}
	}
	
	
	private ArrayList<int[]> getProfilesAtPeaks(List<Region> rs, boolean useCache, String genpath){
		ArrayList<int[]> ret = new ArrayList<int[]>();
		
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genpath);
		}
		
		for(Region r: rs){
			int numK = 0;
			for(int k=minK; k<=maxK; k++){
				numK += (int)Math.pow(4, k);
			}
			int[] pfl = new int[numK];
			String seq = seqgen.execute(r).toUpperCase();
			for(int k=minK; k<=maxK; k++){
				for(int i=0; i<(seq.length()-k+1); i++){
					String currk = seq.substring(i, i+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int currKInt = RegionFileUtilities.seq2int(currk);
					int revcurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revcurrKInt ? currKInt : revcurrKInt;
					int baseind = this.getKmerBaseInd(currk);	
					pfl[baseind+kmer]++;
				}
			}
			ret.add(pfl);
		}
		
		return ret;
	}
	
	
	
	/**
	 * Finds mountains for a given list of regions and given scoring threshold
	 * Odds threshold can be -ve when scanning the model at negPeaks. +ve when scanning the model at posPeaks
	 * Also populates the 
	 * @param useCache
	 * @param genpath
	 * @param oddsThresh
	 */
	private List<Pair<Region,Double>> findMountains(boolean useCache, String genpath, double oddsThresh, List<Region> rs ){
		
		// Detects if scanning is done for the positive or neg set from the sign of oddsThresh
		int classDetector = oddsThresh>0 ? 1:-1;
		
		
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genpath);
		}
		
		List<Pair<Region,Double>> ret = new ArrayList<Pair<Region,Double>>();
		
		for(Region r : rs){
			String seq = seqgen.execute(r).toUpperCase();
			if(seq.contains("N"))
				continue;
			List<Pair<Region,Double>> mountains = new ArrayList<Pair<Region,Double>>();
			for(int l=minM; l<=maxM; l++){
				for(int i=0; i<(seq.length()-l+1); i++){
					String motif = seq.substring(i, i+l);
					double score=0.0;
					for(int k=minK; k<=maxK; k++){
						for(int j=0; j<motif.length()-k+1; j++){
							String currk = motif.substring(j, j+k);
							String revcurrk = SequenceUtils.reverseComplement(currk);
							int  currKInt = RegionFileUtilities.seq2int(currk);
							int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
							int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
							int baseInd = this.getKmerBaseInd(currk);
							score = score+kmerweights[baseInd+kmer];
						}
					}
					if(isPairModel){
						for(int j=0; j<(motif.length()-minK+1-1); j++){
							String currKP1 = motif.substring(j, j+minK);
							String revCurrKP1 = SequenceUtils.reverseComplement(currKP1);
							int currKP1Int = RegionFileUtilities.seq2int(currKP1);
							int revCurrKP1Int = RegionFileUtilities.seq2int(revCurrKP1);
							int kmerP1 = currKP1Int < revCurrKP1Int? currKP1Int : revCurrKP1Int;
							for(int s=j+1;s<(motif.length()-minK+1);s++){
								String currKP2 = motif.substring(s, s+minK);
								String revCurrKP2 = SequenceUtils.reverseComplement(currKP2);
								int currKP2Int = RegionFileUtilities.seq2int(currKP2);
								int revCurrKP2Int = RegionFileUtilities.seq2int(revCurrKP2);
								int kmerP2 = currKP2Int < revCurrKP2Int ? currKP2Int: revCurrKP2Int;
								int x = kmerP1 < kmerP2 ? kmerP1: kmerP2;
								int y = kmerP1 < kmerP2 ? kmerP2: kmerP1;
								score = score+kmerpairwights[x][y];
							}
							
						}
					}
					
					Region hill = new Region(gcon.getGenome(),r.getChrom(),r.getStart()+i,r.getStart()+i+l-1);
					
					Iterator<Pair<Region,Double>> it = mountains.iterator();
					boolean add=true;
					while(it.hasNext() && add){
						Pair<Region,Double> pr = it.next();
						Region currHill = pr.car();
						Double currScore = pr.cdr();
						if(currHill.overlaps(hill) && currScore*classDetector<score*classDetector){
							it.remove();
							add=true;
						}else if(currHill.overlaps(hill) && currScore*classDetector> score*classDetector){
							add=false;
						}
					}
					if(add && score*classDetector > oddsThresh*classDetector){
						mountains.add(new Pair<Region,Double>(hill,score));
					}
					
				}
			}
			ret.addAll(mountains);	
		}
		return ret;
	}
	
	
	
	
	public static void main(String[] args) throws IOException{
		GenomeConfig gcon = new GenomeConfig(args);
		KmerModelScannerCopy scanner = new KmerModelScannerCopy(gcon);
		ArgParser ap = new ArgParser(args);
		
		int mink = Args.parseInteger(args, "minK", 4);
		scanner.setMinK(mink);
		int maxk = Args.parseInteger(args, "maxK", 7);
		scanner.setMaxK(maxk);
		int m = Args.parseInteger(args, "m", 5);
		scanner.setminM(m);
		int M = Args.parseInteger(args, "M", 10);
		scanner.setmaxM(M);
		int win = Args.parseInteger(args, "win", 150);
		
		//String posPeaksFile = ap.getKeyValue("posPeaks");
		
		// Check the string arguments for all posPeaks (THERE CAN BE MULTIPLE posPeaks options)....
		Collection<String> posPeaksFiles = Args.parseStrings(args, "posPeaks");
		List<Point> posPs = new ArrayList<Point>();
		List<Region> posRs = new ArrayList<Region>();
		
		for(String posPeaksFile : posPeaksFiles){
			posPs.addAll(RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), posPeaksFile, win));
			posRs.addAll(RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), posPeaksFile, win));
		}
		
		
		if(ap.hasKey("negPeaks")){
			String negPeaksFile = ap.getKeyValue("negPeaks");
			List<Point> negPs = RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), negPeaksFile, win);
			List<Region> negRs = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), negPeaksFile, win);
			scanner.setNegPeaks(negPs);
		    scanner.setNegRegions(negRs);
		}
			
		int numK = 0;
		for(int k=mink; k<=maxk; k++){
			numK += (int)Math.pow(4, k);
		}
		double[] ws = new double[numK];
		String weightsFile = Args.parseString(args, "weights", "");
		
		double[][] pws=null;
		boolean isPair=false;
		if(ap.hasKey("PairModel")){
			isPair = true;
		}
		BufferedReader reader = new BufferedReader(new FileReader(weightsFile));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\t");
            if(!isPair || !words[0].contains("-")){
            	int kind = scanner.getKmerBaseInd(words[0]) + RegionFileUtilities.seq2int(words[0]);
            	ws[kind] = Double.parseDouble(words[1]);
            }else if(isPair && words[0].contains("-")){
            	String[] subwords = words[0].split("-");
            	int kp1ind = RegionFileUtilities.seq2int(subwords[0]);
            	int kp2ind = RegionFileUtilities.seq2int(subwords[1]);
            	pws = new double[(int) Math.pow(4, mink)][(int) Math.pow(4, mink)];
            	pws[kp1ind][kp2ind] =  Double.parseDouble(words[1]);
            	pws[kp2ind][kp1ind] = Double.parseDouble(words[1]);
            }
        }
        reader.close();
        boolean cache=false;
        String genPath = "";
        if(ap.hasKey("seq")){
        	cache = true;
        	genPath = ap.getKeyValue("seq");
        }
        
        scanner.setKmerWeights(ws);
        if(isPair)
        	scanner.setKmerPairWeights(pws);
        
        String scanType = Args.parseString(args, "scanType", "fixedwin");
        if(!scanType.equals("fixedwin") && !scanType.equals("sliding")){
        	System.err.println("Unrecongnized scanType option; setting to fixedwin type");
        	scanType="fixedwin";
        }
        	
        		
        
        scanner.setModelType(isPair);
        scanner.setPosPeaks(posPs);
        scanner.setPosRegions(posRs);
       
        scanner.setScanType(scanType);
        
        int numClus = Args.parseInteger(args, "numClusters", 3);
        scanner.setNumClusters(numClus);
        
        int kmeansItrs =  Args.parseInteger(args, "kmeansItrs", 100);
        scanner.setClusterItrs(kmeansItrs);
        
        String outbase = Args.parseString(args, "out", "out");
        String wDir = System.getProperty("user.dir");
		String odir = wDir+"/"+outbase;
        File outdir = new File(odir);
        if(outdir.exists())
        	BayesmentsSandbox.deleteDirectory(outdir);
        outdir.mkdir();
        
        scanner.setOutbase(outbase);
        scanner.setOutDir(outdir);
        Double threshold = Double.parseDouble(ap.getKeyValue("oddsthresh"));
        
        scanner.clusterKmerProfilesAtMountains(cache, genPath, threshold);
	}

}
