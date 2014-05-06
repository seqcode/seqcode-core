package edu.psu.compbio.seqcode.projects.shaun.metaexperiment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.motifs.BackgroundModelIO;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

public class DataLoader {
	
	private Genome gen;
	private ArrayList<Region> towers = new ArrayList<Region>();
	//Settings
	private int winSize=200, winStep=100;
	//Points
	private HashMap<Region, Integer> regions = new HashMap<Region, Integer>();
	private HashMap<Integer, Point> points = new HashMap<Integer, Point>();
	private int numRegions=0;
	private boolean entireGenome = false;
	//Data
	private String Yname = new String();
	private double [][] data;
	private String [] names;
	private int numData=0;
	private boolean printBinom=false;
	private boolean needleFiltering=true;
	//Motifs
	private FreqMatrixImport motifImport = new FreqMatrixImport();
	private MarkovBackgroundModel back = null;
	private HashMap<String, WeightMatrix> motifs = new HashMap<String, WeightMatrix>();
	//Experiments
	private HashMap<String, SeqExptInfo> experiments = new HashMap<String, SeqExptInfo>(); 
	
	public static void main(String[] args) throws SQLException, NotFoundException, IOException, ParseException {
		Pair<Organism,Genome> pair = Args.parseGenome(args);
		if(pair==null || !Args.parseArgs(args).contains("y")){printError();return;}
		//if(pair==null){printError();return;}
		
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
		if(Args.parseArgs(args).contains("binomial")){
			loader.setPrintBinomial(true);
		}
		if(Args.parseArgs(args).contains("points")){
		   loader.loadPointsFromFile(Args.parseString(args, "points",null));
		   loader.setEntireGenome(false);
		}if(Args.parseArgs(args).contains("random")){
			loader.loadRandomRegions(Args.parseInteger(args, "random", 10000));
			loader.setEntireGenome(false);
		}
		
		//Data
		//1) Motifs
		if(Args.parseArgs(args).contains("motifs"))
			loader.loadMotifs(Args.parseString(args, "motifs", null), Args.parseString(args, "motifback", null));
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
	public DataLoader(Genome g){
		gen = g;
		System.out.println("DataLoader initialized");
	}
	
	//Accessors
	public void setEntireGenome(boolean e){entireGenome=e;}
	public void setPrintBinomial(boolean e){printBinom=e;}
	public void setWinSize(int w){winSize=w;}
	public boolean examineEntireGenome(){return(entireGenome);}
	public HashMap<String, WeightMatrix> getMotifs(){return motifs;}
	public void setYname(String y){Yname = y;}
	
	
	//Execute
	public void execute(){
		//Check the Y
		if(!motifs.containsKey(Yname) && !experiments.containsKey(Yname)){System.err.println("Y does not exist: "+Yname);}
		
		//Initialize data structures
		for(String s : motifs.keySet()){numData++;}
		for(String s : experiments.keySet()){
			if(experiments.get(s).getType().equals("chip")){numData++;}
		}
		names = new String [numData];
		data = new double[numData][numRegions];
		
		int d=0;
		//Collate the data
		for(String s : motifs.keySet()){
			System.out.println("Motif occupancy for:\t"+s);
			names[d]=s;
			data[d]=calculate(motifs.get(s));
			//Call a scatter-plot generator here
			
			d++;
		}for(String s : experiments.keySet()){
			if(experiments.get(s).getType().equals("chip")){
				System.out.println("ChIP scores for:\t"+s);
				names[d]=s;
				data[d]=calculate(experiments.get(s));
				//Call a scatter-plot generator here
				
				d++;
			}
		}
	}
	
	//Print to a file
	public void printData(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			int y_index=0;	
			//Labels
			fout.write("Region\t"+Yname+"\t");
			for(int x=0; x<names.length; x++){
				if(!names[x].equals(Yname))
					fout.write(names[x]+"\t");
				else
					y_index=x;
			}fout.write("\n");

			//Data
			for(Region r : regions.keySet()){
				int index = regions.get(r);
				fout.write(r.getLocationString()+"\t");
				fout.write(data[y_index][index]+"\t");
				for(int x=0; x<numData; x++){
					if(!names[x].equals(Yname))
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
	//Motif occupancy
	public double [] calculate(WeightMatrix wm){
		double [] scores = new double[numRegions]; for(int i=0; i<numRegions; i++){scores[i]=-1111;}
		WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
		SequenceGenerator seqgen = new SequenceGenerator();
		double conc = Math.exp(-1*wm.getMaxScore());
		if(entireGenome){
			ChromRegionIterator chroms = new ChromRegionIterator(gen);
			while(chroms.hasNext()){
				NamedRegion currentChrom = chroms.next();
				for(int x=currentChrom.getStart(); x<=currentChrom.getEnd(); x+=winStep){
					int y = x+winSize; 
					if(y>currentChrom.getEnd()){y=currentChrom.getEnd();}
					Region r = new Region(gen, currentChrom.getChrom(), x, y);
					String seq = seqgen.execute(r);
					WeightMatrixScoreProfile profiler = scorer.execute(seq);
					double pOcc = profiler2ProbOcc(profiler, conc);
					//scores.add(pOcc);
					//Something...
				}
			}			
		}else{
			for(Region r : regions.keySet()){
				int index = regions.get(r);
				String seq = seqgen.execute(r);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				double pOcc = profiler2ProbOcc(profiler, conc);
				scores[index]=pOcc;
			}
		}
		return(scores);
	}
	//ChIP scores (Z-score if ctrl is defined, read count otherwise)
	public double[] calculate(SeqExptInfo ex){
		double [] scores = new double[numRegions]; for(int i=0; i<numRegions; i++){scores[i]=-1111;}
		SeqExptInfo exp = ex;
		exp.initialize();
		SeqExptInfo ctrl =null;
		if(experiments.containsKey(exp.getCtrl())){
			ctrl = experiments.get(exp.getCtrl());
			ctrl.initialize();
		}
		double ipTotHits = exp.getTotalReads();
		double backTotHits = ctrl==null ? 1 : ctrl.getTotalReads();

		int needleMax = getPoissonThreshold(Math.pow(10, -9), ipTotHits, 1, gen.getGenomeLength(), 0.8, 1, 1);
		
		//Iterate through chromosomes
		ChromRegionIterator chroms = new ChromRegionIterator(gen);
		while(chroms.hasNext()){
			NamedRegion currentChrom = chroms.next();
			//Split the job up into chunks of 50Mbp
			for(int x=currentChrom.getStart(); x<=currentChrom.getEnd(); x+=100000000){
				int y = x+100000000; 
				if(y>currentChrom.getEnd()){y=currentChrom.getEnd();}
				Region currSubRegion = new Region(gen, currentChrom.getChrom(), x, y);
				
				LinkedList<StrandedRegion> ipHits = new LinkedList<StrandedRegion>();
				LinkedList<StrandedRegion> backHits = new LinkedList<StrandedRegion>();
			
				DeepSeqExpt ip = exp.getDeepSeqHandle();
				ipHits.addAll(ip.loadHits(currSubRegion));
				if(ctrl!=null){
					DeepSeqExpt bk = exp.getDeepSeqHandle();
					backHits.addAll(bk.loadHits(currSubRegion));
				}
				if(entireGenome){
					//Do summit
					//Print or store?
				}else{
					for(Region r : regions.keySet()){
						int index = regions.get(r);
						if(r.overlaps(currSubRegion)){
							double ipHitCount=0, backHitCount=0;
							
							ipHitCount = countReadsInWindow(ipHits, r, needleMax, towers);
							if(ctrl!=null && printBinom){
								backHitCount = countReadsInWindow(backHits, r, needleMax, towers);
								scores[index]=binomialSampleEquality(ipHitCount, backHitCount, ipTotHits, backTotHits);
							}else{
								scores[index]=ipHitCount;
							}
						}
					}
				}
			}
		}
		return(scores);
	}
	
	////////////////////////Loaders/////////////////////////////
	//Motifs from freq matrix file & background
	public void loadMotifs(String motifFile, String backFile) throws IOException, ParseException {
		if(motifFile==null){printError(); return;}
		if(backFile==null){
		    back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gen));
        }else{
        	back  = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gen);        	
        }
		motifImport.setBackground(back);
		List<WeightMatrix> motifList = motifImport.readTransfacMatrices(motifFile);
		for(WeightMatrix wm : motifList){
			motifs.put(wm.name, wm);
		}
		//System.out.println(motifs.size());
	}
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
							//System.out.println(potential.getChrom()+":"+potential.getStart()+"-"+potential.getEnd());
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
	            if(words.length>=7){
	            	String name = words[0];
	            	String type = words[1];
	            	String ctrl = words[2];
	            	double rL = new Double(words[3]).intValue();
	            	double rS = new Double(words[4]).intValue();
	            	double z = new Double(words[5]).doubleValue();
	            	SeqExptInfo exp=null;
	            	if(exptType.equals("db")){
		            	ArrayList<SeqLocator> locs = new ArrayList<SeqLocator>();
		            	for(int e=6; e<words.length; e++){
	        		        String[] pieces = words[e].split(";");
	                        if (pieces.length == 2) {
	                            locs.add(new SeqLocator(pieces[0], pieces[1]));
	                        } else if (pieces.length == 3) {
	                            locs.add(new SeqLocator(pieces[0], pieces[1], pieces[2]));
	                        } else {
	                            System.err.println("Couldn't parse a ChipSeqLocator from " + words[e]);
	                        }	                    
		            	}
		            	exp = new SeqExptInfo(gen, name, type, ctrl, rL, rS, locs);
	            	}else{
	            		ArrayList<File> files = new ArrayList<File>();
		            	for(int e=6; e<words.length; e++){
	        		        files.add(new File(words[e]));	                    
		            	}
		            	exp = new SeqExptInfo(gen, name, type, ctrl, rL, rS, files, false, format);
	            	}
	            	exp.setZThres(z);
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
	//Convert a WeightMatrix score to a probability
	private double score2Prob(double x, double conc){
		double kd = Math.exp(-1*x);
		return(conc/(kd+conc));
	}
	//Convert a WeightMatrixProfile into a probability of occupancy 
	private double profiler2ProbOcc(WeightMatrixScoreProfile profile, double conc){
		double prod = 1;
		for(int i=0; i<profile.length(); i++){
			double currScore = profile.getMaxScore(i);
			if(!Double.isNaN(currScore) && !Double.isInfinite(currScore))
				prod*=(1-score2Prob(currScore, conc));
		}
		return(1-prod);
	}
	/* Binomial test for differences between two population proportions */
	protected double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		if(P1==0 && P2==0){return(0);}
		double P = (X1+X2)/(n1+n2);
		double Z0 = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		return(Z0);
	}
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
	private double countReadsInWindow(LinkedList<StrandedRegion> hits, Region reg, int perBaseMax, ArrayList<Region> towers){
		ArrayList<Region> currTowers = new ArrayList<Region>();
		for(Region t : towers)
			if(reg.overlaps(t))
				currTowers.add(t);
		int [] counts = new int[reg.getWidth()+1];
		for(int i=0; i<=reg.getWidth(); i++){counts[i]=0;}
		int count=0;
		for(StrandedRegion r : hits){
			if(r.overlaps(reg)){
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
		}
		return(count);
	}
	
	////////////////////////////////////////////////////////////////////////////
	public static void printError(){
		System.err.println("Usage:\n " +
                "DataLoader \n " +
                "Required:\n  " +
                "--species <organism name;genome version> \n  " +
                "--y <designate something to y>\n " +
                "Options: \n  " +
                "--entiregenome [make regions from whole genome> " +
                "--points <points file> --random <num random points>\n  " +
                "--motifs <motiffile> --motifback <background model>\n  " +
                "--dbexptfile <file listing db experiments to load>\n" +
                "  OR\n" +
                "--exptfile <file listing experiments to load> --format <format of files> ***NO OPTION FOR NON-UNIQUE YET***\n" +
                "--binomial <print binomial test scores where applicable>\n  " +
                "--towers <file with towers to screen out>\n  " +
                "--win <win size>\n  " +
                "--out <output file name>\n  "
                
                );
	}
}
