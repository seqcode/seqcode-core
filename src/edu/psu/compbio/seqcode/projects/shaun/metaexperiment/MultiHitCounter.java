package edu.psu.compbio.seqcode.projects.shaun.metaexperiment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.BackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

public class MultiHitCounter {
	
	private Genome gen;
	private ArrayList<Region> towers = new ArrayList<Region>();
	//Settings
	private int readLength=32;
	private int winSize=500;
	//Points
	private HashMap<Region, Integer> regions = new HashMap<Region, Integer>();
	private HashMap<Integer, Point> points = new HashMap<Integer, Point>();
	private int numRegions=0;
	//Data
	private double [][] data;
	private String [] names;
	private int numData=0;
	private boolean needleFiltering=true;
	//Experiments
	private HashMap<String, DeepSeqExpt> experiments = new HashMap<String, DeepSeqExpt>(); 
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		Pair<Organism,Genome> pair = Args.parseGenome(args);
		if(pair==null){printError();return;}
		
		//Initialize
		DataLoader loader = new DataLoader(pair.cdr());
		
		//Options
		loader.setYname(Args.parseString(args, "y", null));
		String outName = Args.parseString(args, "out", "out.data");
		String towerFileName = Args.parseString(args, "towers", null);
		
		if(Args.parseArgs(args).contains("win")){
			loader.setWinSize(Args.parseInteger(args, "win", 200));
		}
		//Points to load: entire genome, points, or random 
		loader.setEntireGenome(Args.parseFlags(args).contains("entiregenome"));
		//More options
		if(Args.parseArgs(args).contains("points")){
		   loader.loadPointsFromFile(Args.parseString(args, "points",null));
		   loader.setEntireGenome(false);
		}if(Args.parseArgs(args).contains("random")){
			loader.loadRandomRegions(Args.parseInteger(args, "random", 10000));
			loader.setEntireGenome(false);
		}
		
		//Data
		//2) Experiments
		if(Args.parseArgs(args).contains("dbexptfile"))
			loader.loadExperiments(Args.parseString(args, "dbexptfile", null), "db", null);
		if(Args.parseArgs(args).contains("exptfile"))
			loader.loadExperiments(Args.parseString(args, "exptfile", null),"file", Args.parseString(args, "format", "ELAND"));
		
		//Load towers
		if(towerFileName!=null)
			loader.loadTowers(towerFileName);
		
		//Generate the data
		loader.execute();
		
		//Print
		loader.printData(outName);
	}
	
	//Constructor
	public MultiHitCounter(Genome g){
		gen = g;
		System.out.println("DataLoader initialized");
	}
	
	//Accessors
	public void setWinSize(int w){winSize=w;}
	
	
	//Execute
	public void execute(){
		
		//Initialize data structures
		numData = experiments.size();
		names = new String [numData];
		data = new double[numData][numRegions];
		
		int d=0;
		//Collate the data
		for(String s : experiments.keySet()){
			System.out.println("ChIP scores for:\t"+s);
			names[d]=s;
			data[d]=calculate(experiments.get(s));
			
			d++;			
		}
	}
	
	//Print to a file
	public void printData(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			
			//Data
			for(Region r : regions.keySet()){
				int index = regions.get(r);
				fout.write(r.getLocationString()+"\t");
				for(int x=0; x<numData; x++){
					fout.write(data[x][index]+"\t");
				}fout.write("\n");
			}fout.write("\n");
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}
	//Load towers from a file
	public void loadTowers(String tFile){
		towers = loadRegionsFromFile(tFile);
	}
	
	////////////////////////Data Calculators////////////////////
	//ChIP scores (Z-score if ctrl is defined, read count otherwise)
	public double[] calculate(DeepSeqExpt ex){
		double [] scores = new double[numRegions]; for(int i=0; i<numRegions; i++){scores[i]=-1111;}
		DeepSeqExpt exp = ex;
		double ipTotHits = exp.getHitCount();
		
		int needleMax = getPoissonThreshold(Math.pow(10, -9), ipTotHits, 1, gen.getGenomeLength(), 0.8, 1, 1);
		
		for(Region r : regions.keySet()){
			int index = regions.get(r);
			double ipHitCount=0;
			ipHitCount = countReadsInWindow(ex, r, needleMax, towers);
			scores[index]=ipHitCount;
		}
		return(scores);
	}
	
	////////////////////////Loaders/////////////////////////////
	//Points from a file
	public void loadPointsFromFile(String f){
		try {
			File pFile = new File(f);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
			
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            PointParser pparser = new PointParser(gen);
	            Point p=null;
	            if(words.length>=3)
		           	p = pparser.execute(words[2]);
	            else if(words.length==1)
		           	p = pparser.execute(words[0]);
	         	int rstart = p.getLocation()-(winSize/2)<1 ? 1:p.getLocation()-(winSize/2);
	            int rend = p.getLocation()+(winSize/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(winSize/2)-1;
	            Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
	            regions.put(r, numRegions);
	            points.put(numRegions, p);
	            numRegions++;
	        }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Regions from a file
	public ArrayList<Region> loadRegionsFromFile(String f){
		ArrayList<Region> res = new ArrayList<Region>();
		try {
			File pFile = new File(f);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
			
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            RegionParser rparser = new RegionParser(gen);
	            Region r=null;
	            if(words.length>=1)
		           	r = rparser.execute(words[0]);
	         	res.add(r);	            
	        }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(res);
	}
	//Points from random
	public void loadRandomRegions(int numSamples){
		int validSamples=0;
		int genomeSize=0, numChroms=0;
		long [] chromoSize = new long[gen.getChromList().size()];
		Random rand = new Random();
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
					if(randPos+winSize<total+chromoSize[c]){
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+winSize-total));
						
						//is this overlapping an already sampled region?
						boolean over=false;
						for(Region r : regions.keySet()){
							if(r.overlaps(potential)){
								over=true;
							}
						}
						if(!over){
							validSamples++;
							regions.put(potential, numRegions);
							numRegions++;
						}
					}
				}total+=chromoSize[c];
			}				
		}			
	}
	//Load experiment information from a file
	public void loadExperiments(String exptFileName, String exptType, String format){
		try {
			File eFile = new File(exptFileName);
			if(!eFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(eFile));
			String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\t");
	            
	            //read a new SeqExptInfo
	            if(words.length>=2){
	            	String name = words[0];
	            	DeepSeqExpt exp=null;
	            	if(exptType.equals("db")){
		            	ArrayList<ChipSeqLocator> locs = new ArrayList<ChipSeqLocator>();
		            	for(int e=1; e<words.length; e++){
	        		        String[] pieces = words[e].split(";");
	                        if (pieces.length == 2) {
	                            locs.add(new ChipSeqLocator(pieces[0], pieces[1]));
	                        } else if (pieces.length == 3) {
	                            locs.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
	                        } else {
	                            System.err.println("Couldn't parse a ChipSeqLocator from " + words[e]);
	                        }	                    
		            	}
		            	exp = new DeepSeqExpt(gen, locs, "db", readLength);
	            	}else{
	            		ArrayList<File> files = new ArrayList<File>();
		            	for(int e=6; e<words.length; e++){
	        		        files.add(new File(words[e]));	                    
		            	}
		            	exp = new DeepSeqExpt(gen, files, false, format, readLength);
	            	}
	            	experiments.put(name, exp);
	            }            
	        }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/////////////////////////////Extra methods///////////////////////////////////
	//Poisson threshold for needle filter
	public int getPoissonThreshold(double threshold, double count, double hitLength, double seqLen, double mappable, double binWidth, double binOffset){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double pMean = (count*(hitLength/binOffset + binWidth/binOffset))/(seqLen*mappable/binOffset); 
		P.setMean(pMean);
		double l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		return(Math.max(1,countThres));
	}
	// Count the reads in a list that overlap a given region
	// Needle & tower filters
	private double countReadsInWindow(DeepSeqExpt ex, Region reg, int perBaseMax, ArrayList<Region> towers){
		List<ReadHit> hits = ex.loadHits(reg);
		ArrayList<Region> currTowers = new ArrayList<Region>();
		for(Region t : towers)
			if(reg.overlaps(t))
				currTowers.add(t);
		int [] counts = new int[reg.getWidth()+1];
		for(int i=0; i<=reg.getWidth(); i++){counts[i]=0;}
		int count=0;
		for(ReadHit r : hits){
			boolean inTower=false;
			for(Region t : currTowers)
				if(r.overlaps(t))
					inTower=true;
			if(!inTower){
				int offset=r.getStart()-reg.getStart()<0 ? 0 : r.getStart()-reg.getStart();
				if(offset>reg.getWidth())
					offset=reg.getWidth();
				counts[offset]++;
				if(!needleFiltering || (counts[offset] <= perBaseMax)){
					count++;				
				}
			}			
		}
		return(count);
	}
	
	////////////////////////////////////////////////////////////////////////////
	public static void printError(){
		System.err.println("Usage:\n " +
                "MultiHitCounter \n " +
                "Required:\n  " +
                "--species <organism name;genome version> \n  " +
                "Options: \n  " +
                "--points <points file> --random <num random points>\n  " +
                "--dbexptfile <file listing db experiments to load>\n" +
                "  OR\n" +
                "--exptfile <file listing experiments to load> --format <format of files> ***NO OPTION FOR NON-UNIQUE YET***\n" +
                "--towers <file with towers to screen out>\n  " +
                "--win <win size>\n  " +
                "--out <output file name>\n  "
                );
	}
}
