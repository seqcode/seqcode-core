package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.util.Collections;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import cern.jet.stat.Probability;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class MotifAnalysisMultiMotif {

	private GenomeConfig gcon=null;
	private Genome gen =null;
	private List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	private HashMap<String,Double> motifThresholds = new HashMap<String,Double>();
	private MarkovBackgroundModel back;
	private MarkovBackgroundModel simback;
	private List<Region> posSet;
	private List<Point> posPeaks;
	private List<Region> negSet;
	private List<String> posSeq;
	private List<String> negSeq;
	private List<String> posLines;
	private int numRand=1000000;
	private double thresLevel=0.05;
	private double defaultThres=0.0;
	private int window=200;
	private int backOrder=3;
	private boolean printROCCurve=false;
	private int rocStep=100;
	private int rocSlopeWin=1000;
	private int rocSlopeStep=100;
	
	public static void main(String[] args) throws IOException, ParseException {
		ArgParser ap = new ArgParser(args);
		GenomeConfig gConfig = new GenomeConfig(args);
        if(!ap.hasKey("motiffile")||!ap.hasKey("back")) { 
            System.err.println("Usage:\n " +
                               "MotifAnalysisMultiMotif \n" +
                               " Required: \n" +
                               "  --species <organism;genome> " +
                               "  --motiffile <file containing motifs> \n"+
                               "  --back <background Markov model> \n" +
                               " More Information: \n" +
                               "  --seq <path to genome FASTA files>\n" +
                               "  --peaks <file containing coordinates of peaks> \n" +
                               "  --neg <random/markov/filename> \n"+
                               "  --motifthres <file with thresholds> \n" +
                               "  --threslevel <threshold level> OR --multithres \n" +
                               "  --globalthres <fixed threshold for all motifs> \n" +
                               "  --fracthres <fraction of maximum score (all motifs)> \n" +
                               "  --win <window of sequence around positive/negative points> \n"+
                               "  --numrand <number of random sequences to sample> \n" +
                               "  --simback <Markov back for simulating seq>\n" +
                               " Options: \n" +
                               "  --peakswithmotifs [peaks containing ANY motifs] \n" +
                               "  --peaksandmotifs [peaks containing ANY motifs] \n" +
                               "  --peaksandmotifsbest [peaks containing ANY motifs -- best motif only printed] \n" +
                               "  --peakwinmaxscores [maximum motif scores for each peak window] \n" +
                               "  --hitstats [freq/overrep statistics for each motif] \n" +
                               "  --bitpattern [print present/absent bit patterns] \n" +
                               "  --countpattern [print motif count patterns] \n" +
                               "  --printroc" +
                               "");
            return;
        }
    	String motifFile = ap.getKeyValue("motiffile");
    	double globalThreshold =ap.hasKey("globalthres") ? new Double(ap.getKeyValue("globalthres")).doubleValue():Double.MIN_VALUE;  
    	double fractionThreshold =ap.hasKey("fracthres") ? new Double(ap.getKeyValue("fracthres")).doubleValue():Double.MIN_VALUE;
    	String thresFile = ap.hasKey("motifthres") ? ap.getKeyValue("motifthres"):null;
    	double thresLevel = ap.hasKey("threslevel") ? new Double(ap.getKeyValue("threslevel")).doubleValue():0.05;
        String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;
        String simBackFile =ap.hasKey("simback") ? ap.getKeyValue("simback"):backFile;
        String posFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks"):null;
        String neg = ap.hasKey("neg") ? ap.getKeyValue("neg"):null;
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        int numSamp = 1000000;
        if(ap.hasKey("numrand")){
        	numSamp = new Integer(ap.getKeyValue("numrand")).intValue();
        }
        //options
        boolean peaksWithMotifs = ap.hasKey("peakswithmotifs");
        boolean peaksAndMotifs = ap.hasKey("peaksandmotifs");
        boolean peaksAndMotifsBest = ap.hasKey("peaksandmotifsbest");
        boolean peakWinMaxScores = ap.hasKey("peakwinmaxscores");
        boolean hitStats = ap.hasKey("hitStats");
        boolean multiThres = ap.hasKey("multithres");
        boolean bitPattern = ap.hasKey("bitpattern");
        boolean countPattern = ap.hasKey("countpattern");
        boolean print_roc = ap.hasKey("printroc");


		//initialize
		MotifAnalysisMultiMotif analyzer = new MotifAnalysisMultiMotif(gConfig);
		
		//load options
		analyzer.setNumTest(numSamp);
		analyzer.setWin(win);
		analyzer.setPrintROC(print_roc);
		analyzer.loadBackgroundFromFile(backFile, simBackFile);
		analyzer.loadMotifsFromFile(motifFile);
		if(thresFile!=null)
			analyzer.loadThresholdsFromFile(thresFile, thresLevel);
		else if(fractionThreshold != Double.MIN_VALUE)
			analyzer.setAllThresholdsFraction(fractionThreshold);
		else if(globalThreshold != Double.MIN_VALUE)
			analyzer.setAllThresholds(globalThreshold);
					
		//load positive & negative sets
		if(ap.hasKey("seq"))
		{
			analyzer.loadPositive(posFile,true);
			analyzer.loadNegative(neg,true);
		}else{
			analyzer.loadPositive(posFile,false);
			analyzer.loadNegative(neg,false);
		}
		
		

		//Options
		if(peaksWithMotifs)
			analyzer.printPeaksWithMotifs();
		if(peaksAndMotifs)
			analyzer.printBestMotifHits(false);
		if(peaksAndMotifsBest)
			analyzer.printBestMotifHits(true);
		if(peakWinMaxScores)
			analyzer.printPeakWinMaxScores();
		if(hitStats){
			if(multiThres){
				analyzer.loadThresholdsFromFile(thresFile, 0.1);
				analyzer.printHitStats();
				analyzer.loadThresholdsFromFile(thresFile, 0.05);
				analyzer.printHitStats();
				analyzer.loadThresholdsFromFile(thresFile, 0.01);
				analyzer.printHitStats();
				analyzer.loadThresholdsFromFile(thresFile, 0.005);
				analyzer.printHitStats();
				analyzer.loadThresholdsFromFile(thresFile, 0.001);
				analyzer.printHitStats();
			}else{
				analyzer.printHitStats();
			}
		}
		if(bitPattern)
			analyzer.printBitPattern();
		if(countPattern)
			analyzer.printCountPattern();
		
		//analyzer.printMotifInfo();
	}
	
	public MotifAnalysisMultiMotif(GenomeConfig gc){
		gcon = gc;
		gen = gcon.getGenome();
	}
	///////////////////////////////////////////////////////////////////////
	//Options first
	///////////////////////////////////////////////////////////////////////

	//Simple printing of peak lines that contain ANY of the motifs
	public void printPeaksWithMotifs(){
		boolean [] contains = new boolean[posSet.size()];
		for(int i=0; i<posSet.size(); i++){contains[i]=false;}
		for(WeightMatrix m : motifs){
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				if(profiler.getMaxScore()>= motifThresholds.get(m.getName()))
					contains[s]=true;
			}
		}
		for(int i=0; i<posSet.size(); i++){
			if(contains[i])
				System.out.println(posLines.get(i));
		}
	}
	//Print best hits in regions for each motif (only prints if the region contains ANY motif)
	public void printBestMotifHits(boolean printBestOnly){
		boolean [] contains = new boolean[posSet.size()];
		String [][] bestHits = new String[motifs.size()][posSet.size()];
		int [][] bestOffset = new int[motifs.size()][posSet.size()];
		double [][] bestScores = new double[motifs.size()][posSet.size()];
		for(int i=0; i<motifs.size(); i++)
			for(int j=0; j<posSet.size(); j++){
				bestHits[i][j]="NONE";
				bestOffset[i][j]=0;
				bestScores[i][j]=0.0;
			}		
		for(int i=0; i<posSet.size(); i++){contains[i]=false;}
		int x=0;
		for(WeightMatrix m : motifs){
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				if(profiler.getMaxScore()>= motifThresholds.get(m.getName())){
					contains[s]=true;
					bestScores[x][s]=profiler.getMaxScore();
					int index = profiler.getMaxIndex();
					bestOffset[x][s]=Math.abs((posPeaks.get(s).getLocation()-posSet.get(s).getStart())-index);
					String bestSeq = seq.substring(index, index+m.length());
					if(profiler.getMaxStrand()=='-'){
						bestSeq = SequenceUtils.reverseComplement(bestSeq);
					}
					bestHits[x][s]=bestSeq;
				}
			}
			x++;
		}
		for(int i=0; i<posSet.size(); i++){
			if(contains[i]){
				System.out.print(posLines.get(i));
				x=0; int best=0;
				for(WeightMatrix m : motifs){
					if(!printBestOnly)
						System.out.print("\t"+m.getName()+"\t"+bestScores[x][i]+"\t"+bestHits[x][i]);
					if(bestScores[x][i]>bestScores[best][i])
						best=x;
					x++;
				}
				if(printBestOnly)
					System.out.print("\t"+motifs.get(best).getName()+"\t"+bestScores[best][i]+"\t"+bestHits[best][i]);
				System.out.print("\n");
			}
		}
	}
	//Print best hits in regions for each motif
	public void printPeakWinMaxScores(){
		String [][] bestHits = new String[motifs.size()][posSet.size()];
		int [][] bestOffset = new int[motifs.size()][posSet.size()];
		double [][] bestScores = new double[motifs.size()][posSet.size()];
		for(int i=0; i<motifs.size(); i++)
			for(int j=0; j<posSet.size(); j++){
				bestHits[i][j]="NONE";
				bestOffset[i][j]=0;
				bestScores[i][j]=0.0;
			}		
		int x=0;
		for(WeightMatrix m : motifs){
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				bestScores[x][s]=profiler.getMaxScore();
				int index = profiler.getMaxIndex();
				bestOffset[x][s]=Math.abs((posPeaks.get(s).getLocation()-posSet.get(s).getStart())-index);
				String bestSeq = seq.substring(index, index+m.length());
				if(profiler.getMaxStrand()=='-'){
					bestSeq = SequenceUtils.reverseComplement(bestSeq);
				}
				bestHits[x][s]=bestSeq;		
			}
			x++;
		}
		for(int i=0; i<posSet.size(); i++){
			System.out.print(posLines.get(i));
			x=0;
			for(WeightMatrix m : motifs){
				System.out.print("\t"+m.getName()+"\t"+bestScores[x][i]+"\t"+bestHits[x][i]);
				x++;
			}System.out.print("\n");			
		}
	}
	//Print some occurrence and over-representation info for each motif
	public void printHitStats(){
		System.out.println("Threshold: "+thresLevel+"\nMotif\tPosTotal\tPosHits\tPosHitRate\tPosPeaks\tPosPeaksRate\tNegTotal\tNegHits\tNegHitRate\tNegPeaks\tNegPeaksRate\tHitOverRep\tHitsPVal\tPeakOverRep\tPeaksPVal\tROC_AUC");
		ArrayList<MotifStats> stats = new ArrayList<MotifStats>();
		ArrayList<IndexedDouble> hscores = new ArrayList<IndexedDouble>();
		ArrayList<IndexedDouble> pscores = new ArrayList<IndexedDouble>();
		double posTotal = (double)posSeq.size(), negTotal = (double)negSeq.size();
		boolean [] anyHitPos = new boolean [posSeq.size()]; for(int i=0; i<posSeq.size(); i++){anyHitPos[i]=false;}
		boolean [] anyHitNeg = new boolean [negSeq.size()]; for(int i=0; i<negSeq.size(); i++){anyHitNeg[i]=false;}
		int mCount=0;
		for(WeightMatrix m : motifs){
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			ArrayList<Double> posMaxScores = new ArrayList<Double>();
			ArrayList<Double> negMaxScores = new ArrayList<Double>();
			//Counters
			double posHits=0, posPeaks=0;
			double negHits=0, negPeaks=0;
			
			//Positive set
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				posMaxScores.add(profiler.getMaxScore());
				boolean goodPeak =false;
				for(int i=0; i<seq.length(); i++){
					if(profiler.getMaxScore(i)>= motifThresholds.get(m.getName())){
						goodPeak=true;
						posHits++;
					}
				}if(goodPeak){
					posPeaks++;
					anyHitPos[s]=true;
				}
			}
			
			//Negative set
			for(int s=0; s<negSeq.size(); s++){
				String seq = negSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				negMaxScores.add(profiler.getMaxScore());
				boolean goodPeak =false;
				for(int i=0; i<seq.length(); i++){
					if(profiler.getMaxScore(i)>= motifThresholds.get(m.getName())){
						goodPeak=true;
						negHits++;
					}
				}if(goodPeak){
					negPeaks++;
					anyHitNeg[s]=true;
				}
			}
			double roc_auc = calcROCAUC(posMaxScores, negMaxScores, printROCCurve, m);
			
			MotifStats curr = new MotifStats(m.name, posTotal, posHits, posPeaks, negTotal, negHits, negPeaks,roc_auc);
			stats.add(curr);
			//hscores.add(new IndexedDouble(mCount, curr.pvalHits));
			pscores.add(new IndexedDouble(mCount, curr.pvalPeaks));
			mCount++;
		}
		//Collections.sort(hscores);
		Collections.sort(pscores);
		//ArrayList<IndexedDouble> n_hscores = benjaminiHochbergCorrection(hscores);
		ArrayList<IndexedDouble> n_pscores = benjaminiHochbergCorrection(pscores);
		//for(IndexedDouble x : n_hscores){ stats.get(x.id).pvalHits = x.value; }
		for(IndexedDouble x : n_pscores){ stats.get(x.id).pvalPeaks = x.value; }
		
		for(MotifStats m : stats){ m.print(); }

		
		//ANY hits
		double anyPosHit=0, anyNegHit=0;
		for(int i=0; i<posSeq.size(); i++)
			if(anyHitPos[i])
				anyPosHit++;
		for(int i=0; i<negSeq.size(); i++)
			if(anyHitNeg[i])
				anyNegHit++;
		double anyPosRate = anyPosHit/posTotal;
		double anyNegRate = anyNegHit/negTotal;
		double anyOverRep = anyPosRate/anyNegRate;
		System.out.println("ANY\t"+posTotal+"\t\t\t"+anyPosHit+"\t"+anyPosRate+"\t"+negTotal+"\t\t\t"+anyNegHit+"\t"+anyNegRate+"\t\t\t\t"+anyOverRep+"\n");
	}
	
	private double calcROCAUC(ArrayList<Double> posMaxScores, ArrayList<Double> negMaxScores, boolean printROC, WeightMatrix motif) {
		double auc = 0;
		if(posMaxScores.size()==0)
			return 0;
		if(negMaxScores.size()==0)
			return 1;
		ArrayList<LabeledDouble> data = new ArrayList<LabeledDouble>();
		for(Double d : posMaxScores)
			data.add(new LabeledDouble(d, 1));
		for(Double d : negMaxScores)
			data.add(new LabeledDouble(d, 0));
		
		Collections.sort(data);
		double pCount = (double)posMaxScores.size();
		double nCount = (double)negMaxScores.size();
		int x=0;
		double possum=0;
		double lastsn=0;
		double lastfpr=0;
		double lastdval = 10000000;
		if(printROC)
			System.out.println("ROC\t"+motif.getName());
		for(LabeledDouble d : data){
			possum+=d.label;
			if(d.dat!=lastdval){
				double sn = possum/pCount;
				double fp = (x+1)-possum;
				double sp = (nCount-fp)/nCount;
				double fpr=1-sp;
				if(x>0){
						    //Rectangle             //Triangle
					auc += ((fpr-lastfpr)*lastsn) + ((sn-lastsn)*(fpr-lastfpr)/2);
				}
				lastfpr=fpr;
				lastsn = sn;
				if(printROC && x%rocStep==0)
					System.out.println(sn+"\t"+fpr+"\t"+d.dat);
			}
			lastdval = d.dat;
			x++;
		}
		if(printROC){
			//ROC slope analysis
			boolean inflection=false;
			for(int i=(rocSlopeWin/2); i<data.size()-(rocSlopeWin/2) && !inflection; i+=rocSlopeStep){
				double currPos =0;
				for(int j=0; j<rocSlopeWin; j++)
					currPos += data.get(i+j).label;
				if(currPos<rocSlopeWin){
					double slope = (currPos/pCount)/((rocSlopeWin-currPos)/nCount);
					if(slope<1.0){
						inflection=true;
						System.out.println("\n"+motif.getName()+" slope inflection point:\t"+data.get(i).dat);
					}
				}
			}
			if(!inflection)
				System.out.println("\n"+motif.getName()+" No slope inflection point");
		}
		return auc;
	}
	protected class LabeledDouble implements Comparable<LabeledDouble>{
		public Double dat;
		public Integer label;
		public LabeledDouble(Double d, Integer i){dat=d; label=i;}
		public int compareTo(LabeledDouble ld) {
			if(dat > ld.dat){return(-1);}
			else if(dat < ld.dat){return(1);}
			else{return 0;}
		}
	}
	
	//Print a bit patterns for each peak
	public void printBitPattern(){
		boolean [][] contains = new boolean[posSet.size()][motifs.size()];
		for(int i=0; i<posSet.size(); i++)
			for(int j=0; j<motifs.size(); j++)
				contains[i][j]=false;
		 
		for(int j=0; j<motifs.size(); j++){
			WeightMatrix m = motifs.get(j);
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				if(profiler.getMaxScore()>= motifThresholds.get(m.getName()))
					contains[s][j]=true;
			}
		}
		System.out.print("Peak");
		for(int j=0; j<motifs.size(); j++){System.out.print("\t"+motifs.get(j).getName());}
		System.out.print("\n");
		for(int i=0; i<posSet.size(); i++){
			System.out.print(posSet.get(i));
			for(int j=0; j<motifs.size(); j++){
				if(contains[i][j])
					System.out.print("\t1");
				else
					System.out.print("\t0");
			}System.out.print("\n");
		}
	}
	//Print a count patterns for each peak
	public void printCountPattern(){
		int [][] counts = new int[posSet.size()][motifs.size()];
		for(int i=0; i<posSet.size(); i++)
			for(int j=0; j<motifs.size(); j++)
				counts[i][j]=0;
		 
		for(int j=0; j<motifs.size(); j++){
			WeightMatrix m = motifs.get(j);
			WeightMatrixScorer scorer = new WeightMatrixScorer(m);
			
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				WeightMatrixScoreProfile profiler = scorer.execute(seq);
				for(int x=0; x<profiler.length(); x++){
					if(profiler.getMaxScore(x)>= motifThresholds.get(m.getName()))
						counts[s][j]++;	
				}
			}
		}
		System.out.print("Peak");
		for(int j=0; j<motifs.size(); j++){System.out.print("\t"+motifs.get(j).getName());}
		System.out.print("\n");
		for(int i=0; i<posSet.size(); i++){
			System.out.print(posSet.get(i));
			for(int j=0; j<motifs.size(); j++){
				System.out.print("\t"+counts[i][j]);
			}System.out.print("\n");
		}
	}
	///////////////////////////////////////////////////////////////////////
	
	public void setNumTest(int n){numRand=n;}
	public void setWin(int w){window=w;}
	public void setPrintROC(boolean pr){printROCCurve = pr;}

	//load positive
	public void loadPositive(String fname, boolean usecache){
		posSet = RegionFileUtilities.loadRegionsFromPeakFile(gen, fname, window);
		posPeaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, fname, window);
		posLines = RegionFileUtilities.loadLinesFromFile(fname);
		if(usecache){
			posSeq = RegionFileUtilities.getSequencesForRegions(posSet, gcon.getSequenceGenerator());
		}else{
			posSeq = RegionFileUtilities.getSequencesForRegions(posSet, null);
		}
	}
	//load negative
	public void loadNegative(String name, boolean usecache){
		if(name==null || name.equals("random")){
			negSet = RegionFileUtilities.randomRegionPick(gen, posSet, numRand, window);
			negSeq = RegionFileUtilities.getSequencesForRegions(negSet, null);
		}else if(name.equals("markov")){
			negSet = null;
			negSeq = new ArrayList<String>();
			RandomSequenceGenerator rgen = new RandomSequenceGenerator(simback);
			for(int i=0; i<numRand; i++){
				negSeq.add(rgen.execute(window));
			}
		}else{
			negSet = RegionFileUtilities.loadRegionsFromPeakFile(gen, name, window);
			if(usecache){
				negSeq = RegionFileUtilities.getSequencesForRegions(negSet, gcon.getSequenceGenerator());
			}else{
				negSeq = RegionFileUtilities.getSequencesForRegions(negSet, null);
			}
		}
	}
	
	//Load freq matrices
	public void loadMotifsFromFile(String filename){
		FreqMatrixImport motifImport = new FreqMatrixImport();
    	motifImport.setBackground(back);
		motifs.addAll(motifImport.readTransfacMatrices(filename));
		for(WeightMatrix wm : motifs){
			motifThresholds.put(wm.getName(), defaultThres);
		}
	}
	//Load background model
	public void loadBackgroundFromFile(String backFile, String simBackFile) throws IOException, ParseException {		
        if(backFile == null){
        	back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gen));
        }else{
        	back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gen);
        }
        if(!simBackFile.equals(backFile))
        	simback = BackgroundModelIO.parseMarkovBackgroundModel(simBackFile, gen);
        else
        	simback = back;
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
		for(WeightMatrix wm : motifs){
			motifThresholds.put(wm.getName(), t);
		}
	}
	//hard set all thresholds to a fraction of the maximum score
	public void setAllThresholdsFraction(double f){
		for(WeightMatrix wm : motifs){
			motifThresholds.put(wm.getName(), (wm.getMaxScore()*f));
		}
	}
	//Print motif info (testing method)
	public void printMotifInfo(){
		for(WeightMatrix wm : motifs){
			String name = wm.getName();
			int len = wm.length();
			double thres = motifThresholds.get(name);
			
			System.out.println(name+"\t"+len+"\t"+thres);
			System.out.println(wm.printMatrix(wm));
		}
	}
	
	
	//Multiple hypothesis testing correction -- assumes peaks ordered according to p-value
	protected ArrayList<IndexedDouble> benjaminiHochbergCorrection(ArrayList<IndexedDouble> scores){
		double total = scores.size();
		ArrayList<IndexedDouble> res = new ArrayList<IndexedDouble>();
		double rank =1;
		for(IndexedDouble d : scores){
			d.value = d.value*(total/rank);
			if(d.value>1)
				d.value=1.0;
			res.add(new IndexedDouble(d.id, d.value));
			rank++;
		}return(res);
	}
	// Binomial test for differences between two population proportions 
	protected double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		if(!Double.isNaN(Z))
			return(1-Probability.normal(Z));
		else
			return(-1);
	}
	
	protected class MotifStats{
		String name;
		double posTotal, negTotal;
		double posHits=0, posHitRate, posPeaks=0, posPeaksRate;
		double negHits=0, negHitRate, negPeaks=0, negPeaksRate;
		double hitOverRep=0, peaksOverRep=0;
		public double pvalHits=-1, pvalPeaks=-1, ROCAUC = 0.5;
		
		public MotifStats(String name, double posTotal, double posHits, double posPeaks, double negTotal, double negHits, double negPeaks, double rocauc){
			this.name=name; 
			this.posTotal=posTotal;
			this.posHits=posHits; this.posPeaks = posPeaks;
			this.negTotal=negTotal;
			this.negHits=negHits; this.negPeaks = negPeaks;
			this.ROCAUC = rocauc;
			posHitRate = posTotal>0 ? posHits/posTotal : 0;
			posPeaksRate = posTotal>0 ? posPeaks/posTotal : 0;
			negHitRate = negTotal>0 ? negHits/negTotal : 0;
			negPeaksRate = negTotal>0 ? negPeaks/negTotal : 0;
			hitOverRep = negHitRate>0 ? posHitRate/negHitRate:-1;
			peaksOverRep = negPeaksRate>0 ? posPeaksRate/negPeaksRate:-1;
			
			//pvalHits = binomialSampleEquality(posHits, negHits, posTotal, negTotal);
			pvalPeaks = binomialSampleEquality(posPeaks, negPeaks, posTotal, negTotal);
		}

		public void print(){
		    System.out.println(name+"\t"+posTotal+"\t"+posHits+"\t"+posHitRate+"\t"+posPeaks+"\t"+posPeaksRate+"\t"+negTotal+"\t"+negHits+"\t"+negHitRate+"\t"+negPeaks+"\t"+negPeaksRate+"\t"+hitOverRep+"\t"+String.format("%.5e", pvalHits)+"\t"+peaksOverRep+"\t"+String.format("%.5e", pvalPeaks)+"\t"+ROCAUC);
		}
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
