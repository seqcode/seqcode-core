package org.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
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
import org.seqcode.genome.location.ScoredStrandedRegion;
import org.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import org.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gse.gsebricks.verbs.location.PointParser;
import org.seqcode.gse.gsebricks.verbs.location.RegionParser;
import org.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;
import org.seqcode.gse.utils.io.BackgroundModelIO;


public class MultiMotifModuleScanner {

	private Genome gen =null;
	private List<WeightMatrix> primMotifs = new ArrayList<WeightMatrix>();
	private List<WeightMatrix> secMotifs = new ArrayList<WeightMatrix>();
	private HashMap<String,Double> motifThresholds = new HashMap<String,Double>();
	private MarkovBackgroundModel back;
	private List<Region> testRegions;
	private List<Point> testPoints;
	private List<String> testSeq;
	private List<String> testLines;
	private int numRand=1000000;
	private double thresLevel=0.05;
	private double defaultThres=0.0;
	private int window=200;
	
	public static void main(String[] args) throws IOException, ParseException {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("primmotif")||!ap.hasKey("back")) { 
            System.err.println("Usage:\n " +
                               "MotifAnalysisMultiMotif: \n" +
                               "\tPrint peaks containing ANY primary motifs AND ANY secondary motifs\n\n" +
                               " Required: \n" +
                               "  --species <organism;genome> " +
                               "  --primmotif <file containing primary motifs> \n"+
                               "  --back <background Markov model> \n" +
                               " More Information: \n" +
                               "  --secmotif <file containing secondary motifs> \n"+
                               "  --peaks <file containing coordinates of peaks> \n" +
                               "  --motifthres <file with thresholds> \n" +
                               "  --threslevel <threshold level> OR --multithres \n" +
                               "  --globalthres <fixed threshold for all motifs> \n" +
                               "  --fracthres <fraction of maximum score (all motifs)> \n" +
                               "  --win <window of sequence around positive/negative points> \n"+
                               "  --numrand <number of random sequences to sample> \n" +
                               "  --global [get all sliding windows that match the module]\n" +
                               "  --random [scan random windows for modules]\n" +
                               "  --modinfo [print module info instead of scanning]\n" +
                               //" Options: \n" +
                               //"  --peakswithmodules [peaks containing ANY primary motifs AND ANY secondary motifs] \n" +
                               "");
            return;
        }
        try {
			Pair<Species, Genome> pair = Args.parseGenome(args);
			Genome currgen = pair.cdr();
			Collection<String> primMotifFiles = Args.parseStrings(args, "primmotif");
			Collection<String> secMotifFiles = Args.parseStrings(args, "secmotif");
	    	double globalThreshold =ap.hasKey("globalthres") ? new Double(ap.getKeyValue("globalthres")).doubleValue():Double.MIN_VALUE;  
	    	double fractionThreshold =ap.hasKey("fracthres") ? new Double(ap.getKeyValue("fracthres")).doubleValue():Double.MIN_VALUE;
	    	String thresFile = ap.hasKey("motifthres") ? ap.getKeyValue("motifthres"):null;
	    	double thresLevel = ap.hasKey("threslevel") ? new Double(ap.getKeyValue("threslevel")).doubleValue():0.05;
	        String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;
	        String testFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks"):null;
	        //String neg = ap.hasKey("neg") ? ap.getKeyValue("neg"):null;
	        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():200;
	        int numSamp = 1000000;
	        if(ap.hasKey("numrand")){
	        	numSamp = new Integer(ap.getKeyValue("numrand")).intValue();
	        }
	        //options
	       // boolean peaksWithModules = ap.hasKey("peakswithmodules");
	        boolean globalModules = ap.hasKey("global");
	        boolean randomModules = ap.hasKey("random");
	        boolean printModInfo = ap.hasKey("modinfo");
	        

			//initialize
	        MultiMotifModuleScanner analyzer = new MultiMotifModuleScanner(currgen);
			
			//load options
			analyzer.setNumRandom(numSamp);
			analyzer.setWin(win);
			analyzer.loadBackgroundFromFile(backFile);
			if(primMotifFiles!=null)
				for(String s : primMotifFiles)
					analyzer.setPrimMotifs(analyzer.loadMotifsFromFile(s));
			if(secMotifFiles!=null)
				for(String s : secMotifFiles)
					analyzer.setSecMotifs(analyzer.loadMotifsFromFile(s));
			if(thresFile!=null)
				analyzer.loadThresholdsFromFile(thresFile, thresLevel);
			else if(fractionThreshold != Double.MIN_VALUE)
				analyzer.setAllThresholdsFraction(fractionThreshold);
			else if(globalThreshold != Double.MIN_VALUE)
				analyzer.setAllThresholds(globalThreshold);
						
			if(testFile!=null || randomModules){
				if(testFile!=null)
					analyzer.loadTest(testFile);
				if(randomModules)
					analyzer.loadRandom();
				if(printModInfo)
					analyzer.printModuleInfo();
				else
					analyzer.printPeaksWithModules();
			}else if(globalModules){
				analyzer.printModuleRegionsFullGenome();
			}
			
			//analyzer.printMotifInfo();			
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public MultiMotifModuleScanner(Genome g){
		gen=g;
	}
	///////////////////////////////////////////////////////////////////////
	//Options first
	///////////////////////////////////////////////////////////////////////
	public void setPrimMotifs(List<WeightMatrix> w){primMotifs.addAll(w);}
	public void setSecMotifs(List<WeightMatrix> w){secMotifs.addAll(w);}
	public void setNumRandom(int n){numRand=n;}
	public void setWin(int w){window=w;}
	
	//load positive
	public void loadTest(String fname){
		testRegions = loadRegionsFromPeakFile(fname, window);
		testPoints = loadPeaksFromPeakFile(fname, window);
		testLines = loadLinesFromFile(fname);
		testSeq = getSequencesForRegions(testRegions);
	}
	//load random
	public void loadRandom(){
		testRegions = randomRegionPick(null, numRand, window);
		testPoints = regions2midpoints(testRegions);
		testLines = regions2lines(testRegions);
		testSeq = getSequencesForRegions(testRegions);
	}
	
	//Simple printing of peak lines that contain ANY of the primary motifs AND ANY of the secondary motifs
	public void printPeaksWithModules(){
		boolean [] containsPrim = new boolean[testRegions.size()];
		boolean [] containsSec = new boolean[testRegions.size()];
		for(int i=0; i<testRegions.size(); i++){containsPrim[i]=false;}
		for(int i=0; i<testRegions.size(); i++){
			if(secMotifs.size()==0)
				containsSec[i]=true;
			else
				containsSec[i]=false;
		}
		//Primary
		for(WeightMatrix m : primMotifs){
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			for(int s=0; s<testSeq.size(); s++){
				String seq = testSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				if(profiler.getMaxScore()>= motifThresholds.get(m.getName()))
					containsPrim[s]=true;
			}
		}
		//Secondary
		for(WeightMatrix m : secMotifs){
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			for(int s=0; s<testSeq.size(); s++){
				String seq = testSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				if(profiler.getMaxScore()>= motifThresholds.get(m.getName()))
					containsSec[s]=true;
			}
		}
		for(int i=0; i<testRegions.size(); i++){
			if(containsPrim[i] && containsSec[i])
				System.out.println(testLines.get(i));
		}
	}
	
	//printing the hit sequences
	public void printModuleRegionsFullGenome(){
		List<Region> chromRegs = new ArrayList<Region>();
		ChromRegionIterator chroms = new ChromRegionIterator(gen);
		while(chroms.hasNext()){ NamedRegion c = chroms.next(); chromRegs.add(c);}
		SequenceGenerator seqgen = new SequenceGenerator();
		ArrayList<ScoredStrandedRegion> hits = new ArrayList<ScoredStrandedRegion>();
		ArrayList<String> hitseqs = new ArrayList<String>();
		
		for(int i=0; i<chromRegs.size(); i++){
			Region r = chromRegs.get(i);
			
			int rstart = r.getStart();
            int chunksize = 10000000;
            int length = window;
            // work over the target region in pieces
            while (rstart < r.getEnd()) {
                int rend = rstart + chunksize;
                if (rend > r.getEnd()) {
                    rend = r.getEnd();
                }if (rend - rstart < length) {break;}
				
                Region cr = new Region(r.getGenome(), r.getChrom(), rstart, rend);
				String seq = seqgen.execute(cr);
				
				//Find modules in sliding windows
				int winSize = window;
				int winStep = window/2;
				int numWins = (chunksize/winStep)+1;
				boolean [] containsPrim = new boolean[numWins];
				boolean [] containsSec = new boolean[numWins];
				for(int x=0; x<numWins; x++){containsPrim[x]=false;}
				for(int x=0; x<numWins; x++){
					if(secMotifs.size()==0)
						containsSec[x]=true;
					else
						containsSec[x]=false;
				}
			
				//Primary
				for(WeightMatrix m : primMotifs){
					double currThres = motifThresholds.get(m.getName());
					WeightMatrixScorer scorer = new WeightMatrixScorer(m);
					WeightMatrixScoreProfile profiler = scorer.execute(seq);
					for(int s=0; s<numWins; s++){
						int wStart = rstart+(s*winStep);
						int wEnd = wStart + winSize - m.length()+1;
						for(int z=wStart; z<wEnd; z++){
							if(profiler.getMaxScore(z)>= currThres)
								containsPrim[s]=true;
						}
					}
				}
				//Secondary
				for(WeightMatrix m : secMotifs){
					double currThres = motifThresholds.get(m.getName());
					WeightMatrixScorer scorer = new WeightMatrixScorer(m);
					WeightMatrixScoreProfile profiler = scorer.execute(seq);
					for(int s=0; s<numWins; s++){
						int wStart = rstart+(s*winStep);
						int wEnd = wStart + winSize - m.length()+1;
						for(int z=wStart; z<wEnd; z++){
							if(profiler.getMaxScore(z)>= currThres)
								containsSec[s]=true;
						}
					}
				}

				//Print the matching windows
				for(int s=0; s<numWins; s++){
					if(containsPrim[s] && containsSec[s]){
						int wStart = rstart+(s*winStep);
						Region currWin = new Region(gen, r.getChrom(), wStart, wStart+winSize-1);
						System.out.println(currWin);
					}
				}
				
				rstart = rend - length + 1;
            }
		}
	}
	
	//Print module info (motif max scores and number of non-overlapping motif instances over threshold)
	public void printModuleInfo(){
		
		System.out.print("Sequence");
		for(WeightMatrix m : primMotifs)
			System.out.print("\t"+m.name+":count\t"+m.name+":max");
		for(WeightMatrix m : secMotifs)
			System.out.print("\t"+m.name+":count\t"+m.name+":max");
		System.out.print("\n");
		
		for(int s=0; s<testSeq.size(); s++){
			System.out.print(testRegions.get(s));
			String seq = testSeq.get(s);
			//Primary
			for(WeightMatrix m : primMotifs){
				WeightMatrixScorer scorer = new WeightMatrixScorer(m);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				double maxScore = profiler.getMaxScore();
				int count=0;
				for(int x=0; x<seq.length(); x++){
					if(profiler.getMaxScore(x)>=motifThresholds.get(m.getName())){
						count++; 
						x+=m.length();
					}
				}
				System.out.print("\t"+count+"\t"+maxScore);
			}
			//Secondary
			for(WeightMatrix m : secMotifs){
				WeightMatrixScorer scorer = new WeightMatrixScorer(m);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				double maxScore = profiler.getMaxScore();
				int count=0;
				for(int x=0; x<seq.length(); x++){
					if(profiler.getMaxScore(x)>=motifThresholds.get(m.getName())){
						count++;
						x+=m.length();
					}
				}
				System.out.print("\t"+count+"\t"+maxScore);
			}
			System.out.print("\n");
		}
	}
	
	//Get sequences for a set of regions
	private ArrayList<String> getSequencesForRegions(List<Region> regions){
		ArrayList<String> seqs = new ArrayList<String>(); 
		SequenceGenerator seqgen = new SequenceGenerator();
		for(Region r : regions){
			seqs.add(seqgen.execute(r));
		}return(seqs);
	}
	//Load freq matrices
	public List<WeightMatrix> loadMotifsFromFile(String filename){
		FreqMatrixImport motifImport = new FreqMatrixImport();
    	motifImport.setBackground(back);
    	List<WeightMatrix> currmotifs = new ArrayList<WeightMatrix>();
		currmotifs.addAll(motifImport.readTransfacMatrices(filename));
		for(WeightMatrix wm : currmotifs){
			motifThresholds.put(wm.getName(), defaultThres);
		}
		return currmotifs;
	}
	//Load background model
	public void  loadBackgroundFromFile(String backFile) throws IOException, ParseException {		
        if(backFile == null){
        	back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gen));
        }else{
        	back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gen);
        }
	}
	//Load thresholds
	public void loadThresholdsFromFile(String filename, double level){
		thresLevel=level;
		if(filename == null){System.err.println("No threshold file specified");}
		else{
			int thresIndex=3;
			try{
				File bFile = new File(filename);
				if(bFile.isFile()){
					BufferedReader reader;
					reader = new BufferedReader(new FileReader(bFile));
					String firstLine = reader.readLine();
					String [] tokens = firstLine.split("[\\s*\\t\\r\\n\\f]");
					for(int i=3; i<tokens.length; i++){
						String t = tokens[i];
						if(t.startsWith("Thres")){
							double val = new Double(t.replaceAll("Thres", "")).doubleValue();
							if(val==thresLevel){
								thresIndex=i;
							}
						}
					}
					String line;
					while((line= reader.readLine())!=null){
						tokens = line.split("[\\s*\\t\\r\\n\\f]");
						String name = tokens[0];
						double v = new Double(tokens[thresIndex]);
						motifThresholds.put(name, v);
					}
					reader.close();
				}
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}
	}
	//hard set all thresholds
	public void setAllThresholds(double t){
		for(WeightMatrix wm : primMotifs){
			motifThresholds.put(wm.getName(), t);
		}
		for(WeightMatrix wm : secMotifs){
			motifThresholds.put(wm.getName(), t);
		}
	}
	//hard set all thresholds to a fraction of the maximum score
	public void setAllThresholdsFraction(double f){
		for(WeightMatrix wm : primMotifs){
			motifThresholds.put(wm.getName(), (wm.getMaxScore()*f));
		}
		for(WeightMatrix wm : secMotifs){
			motifThresholds.put(wm.getName(), (wm.getMaxScore()*f));
		}
	}
	//Print motif info (testing method)
	public void printMotifInfo(){
		for(WeightMatrix wm : primMotifs){
			String name = wm.getName();
			int len = wm.length();
			double thres = motifThresholds.get(name);
			
			System.out.println(name+"\t"+len+"\t"+thres);
			System.out.println(wm.printMatrix(wm));
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
	            
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		                if(win==-1 && words[0].contains(":") && words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	regs.add(q);
		                }else{
		                	regs.add(p.expand(win/2));
		                }
	                }else if(words.length>=1 && words[0].contains(":")){
		            	if(words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	if(win==-1){
			                	if(q!=null){regs.add(q);}
			                }else
			                	regs.add(q.getMidpoint().expand(win/2));
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	regs.add(p.expand(win/2));
		            	}
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
	//Load a set of regions from a peak file
	public ArrayList<Point> loadPeaksFromPeakFile(String filename, int win){
		ArrayList<Point> peaks = new ArrayList<Point>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	peaks.add(p);		                
	                }else if(words.length>=1 && words[0].contains(":")){
		            	if(words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	peaks.add(q.getMidpoint());			            	
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	peaks.add(p);
		            	}
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
		return(peaks);
	}
	private ArrayList<String> loadLinesFromFile(String filename){
		ArrayList<String> lines = new ArrayList<String>();
		try{
			File inFile = new File(filename);
			if(!inFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(inFile));
	        String line;//Ignore first line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	            	lines.add(line);
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return(lines);
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
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-1-total));
						
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
	
	protected List<Point> regions2midpoints(List<Region> regs){
		List<Point> p = new ArrayList<Point>();
		for(Region r : regs)
			p.add(r.getMidpoint());
		return p;
	}
	
	protected List<String> regions2lines(List<Region> regs){
		List<String> p = new ArrayList<String>();
		for(Region r : regs)
			p.add(r.toString());
		return p;
	}
	
	protected class IndexedDouble implements Comparable<IndexedDouble>{
		public Integer id;
		public Double value;
		
		public IndexedDouble(Integer i, Double v){id=i; value=v;}
		
		public int compareTo(IndexedDouble x) {
		    if(value<x.value){return(-1);}
		    else if(value>x.value){return(1);}
		    else{return(0);}
		}
    }
}
