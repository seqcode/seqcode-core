package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.RepeatMaskedRegion;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.chipseq.MACSParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.chipseq.MACSPeakRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.projects.gps.ReadHit;
import edu.psu.compbio.seqcode.gse.projects.gps.features.Feature;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class PeaksAnalysis {
	private Species org=null;
	private Genome gen =null;
	private List<Region> posSet;
	private List<Point> posPeaks=null;
	private List<String> posLines;
	private List<Region> negSet;
	private List<Point> negPeaks=null;
	private List<String> negLines;
	private List<Region> towers = new ArrayList<Region>();
	private int window=200;
	private boolean shiftTags = false;
	private int readShift=0;
	private int readExtension=0;
	protected String exptLabel;
	protected DeepSeqExpt IP=null;
	protected DeepSeqExpt IP2=null;
	protected boolean dbconnected=false;
	protected boolean twoExpts=false;
	protected boolean normCounts = false;
	private double iphittot=0;
	private int metaBinSize=20;
	private int fixedPBCutoff=-1;
	protected int readLen=32;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") ||(!ap.hasKey("peaks")&&!ap.hasKey("macspeaks")&&!ap.hasKey("motifhits"))) { 
            System.err.println("Usage:\n " +
                               "PeaksAnalysis \n" +
                               " Required: \n" +
                               "  --species <species;version> " +
                               "  --peaks <file containing coordinates of peaks in Shaun's format> \n" +
                               "     OR" +
                               "  --macspeaks <file containing coordinates of peaks in MACS format> \n" +
                               "     OR" +
                               "  --motifhits <file containing motif hit coordinates> \n" +
                               "  --neg <file containing negative regions in Shaun's format> \n" +
                               "     OR" +
                               "  --macsneg <file containing negative regions in MACS format> \n" +
                               "  --towers <file containing towers in Shaun's format>\n" +
                               " More Info: \n"+
                               "  --win <window of sequence around peaks> \n"+
                               " Options: \n" +
                               "  --out output filename\n" +
                               "  --cache [flag to cache sequences] AND --seq <path to genome files>\n" +
                               "  --printseqs [flag to print sequences under peaks] \n" +
                               "  --printseqkmers [flag to print k-mer counts for sequences under peaks] \n" +
                               "  --printxy [flag to print hit counts under peaks for two expts]\n" +
                               "  --readlen <rlen>\n" +
                               "  --fixedpb <cutoff>\n" +
                               "  --expt <IP expt names> \n" +
                               "  --expt2 <IP2 expt names> \n" +
                               "  --exptlabel <experiment name> \n" +
                               "  --format <ELAND/NOVO/BOWTIE>\n" +
                               "  --metapeak [flag to print a metapeak (IP expt handles reqd)] \n" +
                               "  --weightcount <provide a metapeak file, get back weighted counts in the window around peaks>\n" +
                               "  --counts [flag to print raw counts per peak]\n" +
                               "  --orderbycounts [flag to print peaks ordered by raw counts per peak]\n" +
                               "  --quick [flag to do things quickly without tower & needle filters]\n" +
                               "  --normcounts [sum the counts per peak, weighting each experiment by the total number of reads]\n" +
                               "  --repscreen [when printing sequences, screen out anything overlapping a repeat]\n" +
                               "");
            return;
        }
        try {
        	Pair<Species, Genome> pair = Args.parseGenome(args);
        	Species currorg = pair.car();
        	Genome currgen = pair.cdr();
			
	        String peaksFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks") : (ap.hasKey("macspeaks") ? ap.getKeyValue("macspeaks"):null);
	        String motifHitsFile = ap.hasKey("motifhits") ? ap.getKeyValue("motifhits") : null;
	    	boolean macsformatpeaks = ap.hasKey("macspeaks");
	    	String negFile = ap.hasKey("neg") ? ap.getKeyValue("neg") : (ap.hasKey("macsneg") ? ap.getKeyValue("macsneg"):null);
	    	boolean hasNeg = negFile==null ? false : true;
	    	boolean macsformatneg = ap.hasKey("macsneg");
	    	String towerFile = ap.hasKey("towers") ? ap.getKeyValue("towers") : null;
	    	boolean hasTowers = towerFile==null ? false : true;
	    	boolean printxy = ap.hasKey("printxy");
	        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
	        int fpb = ap.hasKey("fixedpb") ? new Integer(ap.getKeyValue("fixedpb")).intValue():-1;
	        int rL = ap.hasKey("readlen") ? new Integer(ap.getKeyValue("readlen")).intValue():32;
	        String outFile = ap.hasKey("out") ? ap.getKeyValue("out") : "out.txt";
	        String label = ap.hasKey("exptlabel") ? ap.getKeyValue("exptlabel") : "IP";
	     	String mFile = ap.hasKey("weightcount") ? ap.getKeyValue("weightcount") : null;
	     	
	        //options
	        boolean printSeqs = ap.hasKey("printseqs");
	        boolean printSeqKmers = ap.hasKey("printseqkmers");
	        boolean printTiledSeqKmers = ap.hasKey("printtiledseqkmers");
	        boolean printSeqKmersCoOccurence = ap.hasKey("printSeqKmersCoOccurence");
	        int printSeqKmersCoOccurenceK = ap.hasKey("printSeqKmersCoOccurence") ? new Integer(ap.getKeyValue("printSeqKmersCoOccurence")).intValue():4;
	        int printSeqKmersK = ap.hasKey("printseqkmers") ? new Integer(ap.getKeyValue("printseqkmers")).intValue():4;;
	        boolean printMeta = ap.hasKey("metapeak");
	        boolean weightedCounts = ap.hasKey("weightcount");
	        boolean normCounts = ap.hasKey("normcounts");
	        boolean counts = ap.hasKey("counts");
	        boolean orderByCounts = ap.hasKey("orderbycounts");
	        boolean quick = ap.hasKey("quick");
	        boolean repeatScreen = ap.hasKey("repscreen");
	        boolean useCache = ap.hasKey("cache");
	        String genomePath = "";
	        if(useCache){
	        	genomePath = ap.getKeyValue("seq");
	        }
	        
        	/////////////////////////////////////////////////////////////////////////
	        ///////////// START 
        
			//initialize
			PeaksAnalysis analyzer = new PeaksAnalysis(currorg, currgen);
			analyzer.setWin(win);
			analyzer.setReadLen(rL);
			analyzer.setFixedPB(fpb);
			
			//load positive & negative sets
			analyzer.loadPeaks(peaksFile, motifHitsFile, macsformatpeaks);
			if(hasNeg)
				analyzer.loadNegs(negFile, macsformatneg);
			if(hasTowers)
				analyzer.loadTowers(towerFile);
			
			//Options
			if(printSeqs)
				analyzer.printPeakSeqs(repeatScreen);
			if(printSeqKmers)
				analyzer.printPeakSeqKmers(printSeqKmersK, useCache, genomePath);
			
			if(printTiledSeqKmers){
				int kmin = ap.hasKey("kmin") ? new Integer(ap.getKeyValue("kmin")).intValue():3;
				int kmax =  ap.hasKey("kmax") ? new Integer(ap.getKeyValue("kmax")).intValue():8;
				analyzer.printPeakSeqTiledKmers(kmin, kmax, useCache, genomePath);
			}
			
			if(printSeqKmersCoOccurence){
				int s = ap.hasKey("slide") ? new Integer(ap.getKeyValue("slide")).intValue():10;
				analyzer.printPeakSeqKmerCoocurrence(printSeqKmersCoOccurenceK, useCache,genomePath , s);
			}
			
			//analyzer.printMotifInfo();
			
			analyzer.loadExperiments(args);
			
			if(printMeta){
				analyzer.buildMetaPeak(outFile);
			}
			if(printxy){
				analyzer.printXY();
			}
			
			if(weightedCounts){
				double [] meta = analyzer.loadMetaPeakFile(mFile);
				analyzer.weightedCounts(meta, label, outFile);
			}
			if(orderByCounts){
				analyzer.orderByCounts(label, outFile);
			}
			if(counts){
				if(normCounts)
					analyzer.setNormCounts(true);
				if(quick)
					analyzer.quickcounts(label, outFile);
				else
					analyzer.counts(label, outFile);
			}
			analyzer.close();
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public PeaksAnalysis(Species o, Genome g){
		org = o;
		gen=g;
	}
	///////////////////////////////////////////////////////////////////////
	//Options first
	///////////////////////////////////////////////////////////////////////

	//Print the sequences for each peak 
	public void printPeakSeqs(boolean repeatScreen){
		RepeatMaskedGenerator repMask = new RepeatMaskedGenerator(gen);
		SequenceGenerator seqgen = new SequenceGenerator();
		for(Region r : posSet){
			String seq = seqgen.execute(r);
			
			if(repeatScreen){
				boolean repOver = false;
				Iterator<RepeatMaskedRegion> repItr = repMask.execute(r);
                while(repItr.hasNext()){
                    RepeatMaskedRegion currRep = repItr.next();
                    if(currRep.overlaps(r)){
                        repOver=true;
                    }
                }
                if(!repOver){
                	if(r instanceof StrandedRegion && ((StrandedRegion)r).getStrand()=='-'){
                		seq = SequenceUtils.reverseComplement(seq);
        	        	System.out.println(">"+r.toString()+"\n"+seq);
                	}else
                		System.out.println(">"+r.toString()+"\n"+seq);
                }
			}else{
				if(r instanceof StrandedRegion && ((StrandedRegion)r).getStrand()=='-'){
            		seq = SequenceUtils.reverseComplement(seq);
    	        	System.out.println(">"+r.toString()+"\n"+seq);
            	}else
            		System.out.println(">"+r.toString()+"\n"+seq);
			}
		}
	}
	
	//Print the k-mers in the sequences for each peak 
	public void printPeakSeqKmers(int k, boolean useCache, String genPath ){
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genPath);
		}
		int numK = (int)Math.pow(4, k);
		int [] kmerCounts = new int[numK];
		System.out.print("Region");
		for(int i=0; i<numK; i++)
			System.out.print("\t"+RegionFileUtilities.int2seq(i, k));
		System.out.println("");
		for(Region r : posSet){
			for(int i=0; i<numK; i++)
				kmerCounts[i]=0;
			
			String seq = seqgen.execute(r).toUpperCase();
			//Check if the sequence (seq) contains any N's if present ignore them
			if(seq.contains("N"))
				continue;
			for(int i=0; i<(seq.length()-k+1); i++){
				String currK = seq.substring(i, i+k);
				String revCurrK =SequenceUtils.reverseComplement(currK);
				int  currKInt = RegionFileUtilities.seq2int(currK);
				int  revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
				int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
				kmerCounts[kmer]++;
			}
			System.out.print(r.getLocationString());
			for(int i=0; i<numK; i++)
				System.out.print("\t"+kmerCounts[i]);
			System.out.println("");
		}
	}
	
	
	// Print the number of times 2 kmers co-occurr in a given sliding window 
	// 									and
	// Prints the kmer frequency counts
	public void printPeakSeqKmerCoocurrence(int k, boolean useCache, String genPath, int s){
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genPath);
		}
		int numK = (int)Math.pow(4, k);
		int[] kmerCounts = new int[numK];
		int[][] kmerCoO = new int[numK][numK];
		
		//Print header
		System.out.print("Region");
		for(int i=0; i<numK; i++){
			System.out.print("\t"+RegionFileUtilities.int2seq(i, k));
		}
		for(int r=0; r<numK; r++){
			for(int c=r; c<numK; c++){
				System.out.print("\t"+RegionFileUtilities.int2seq(r, k)+"-"+RegionFileUtilities.int2seq(c, k));
			}
		}
		System.out.println("");
		
		for(Region r: posSet){
			// Initialize 
			for(int i=0; i<numK; i++)
				kmerCounts[i]=0;
			for(int i=0; i<numK; i++){
				for(int j=0; j<numK; j++){
					kmerCoO[i][j]=0;
				}
			}
			
			String seq = seqgen.execute(r).toUpperCase();
			if(seq.contains("N"))
				continue;
			
			// Count kmer freqs
			for(int i=0; i<(seq.length()-k+1); i++){
				String currK = seq.substring(i, i+k);
				String revCurrK =SequenceUtils.reverseComplement(currK);
				int  currKInt = RegionFileUtilities.seq2int(currK);
				int  revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
				int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
				kmerCounts[kmer]++;
			}
			
			// Count the co-occurrence
			for(int i=0; i<(seq.length()-s+1);i++){
				String currSubSeq = seq.substring(i,i+s);
				for(int j=0;j<(currSubSeq.length()-k+1-1);j++){
					String currKP1 = currSubSeq.substring(j, j+k);
					String revCurrKP1 = SequenceUtils.reverseComplement(currKP1);
					int currKP1Int = RegionFileUtilities.seq2int(currKP1);
					int revCurrKP1Int = RegionFileUtilities.seq2int(revCurrKP1);
					int kmerP1 = currKP1Int<revCurrKP1Int? currKP1Int:revCurrKP1Int;
					for(int l=j+1;l<(currSubSeq.length()-k+1); l++){
						String currKP2 = currSubSeq.substring(l,l+k);
						String revCurrKP2 = SequenceUtils.reverseComplement(currKP2);
						int currKP2Int = RegionFileUtilities.seq2int(currKP2);
						int revCurrKP2Int = RegionFileUtilities.seq2int(revCurrKP2);
						int kmerP2 = currKP2Int<revCurrKP2Int? currKP2Int: revCurrKP2Int;
						kmerCoO[kmerP1][kmerP2]=1;
						kmerCoO[kmerP2][kmerP1]=1;
					}
				}
			}
			System.out.print(r.getLocationString());
			
			for(int i=0; i<numK; i++)
				System.out.print("\t"+kmerCounts[i]);
			for(int i=0; i<numK; i++){
				for(int j=i; j<numK; j++){
					System.out.print("\t"+kmerCoO[i][j]);
				}
			}
			System.out.println("");
		}
	}
	
	
	//Divide the sequence into 3 regions and print kmer counts in the left, right and central regions
	public void printPeakSeqTiledKmers(int kmin, int kmax, boolean useCache, String genPath){
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genPath);
		}
		
		int TnumK = 0;
		for(int k=kmin; k<=kmax; k++){
			TnumK = TnumK + ((int)Math.pow(4, k))*3;
		}
		int[] kmerCounts = new int[TnumK];
		System.out.print("Region");
		int currIndex = 0;
		
		//Print header
		for(int k=kmin; k<=kmax; k++){
			int numK = (int)Math.pow(4, k);
			
			//Print left kmers 
			for(int i=0; i<numK; i++){
				System.out.print("\t"+RegionFileUtilities.int2seq(i, k)+"_L");
			}
			//Print central kmers
			for(int i=0; i<numK; i++){
				System.out.print("\t"+RegionFileUtilities.int2seq(i, k)+"_C");
			}
			// Print right kmers
			for(int i=0; i<numK; i++){
				System.out.print("\t"+RegionFileUtilities.int2seq(i, k)+"_R");
			}
		}
		System.out.println("");
		
		for(Region r : posSet){
			for(int i=0; i<kmerCounts.length; i++){
				kmerCounts[i] = 0;
			}
			
			StrandedRegion sr = (StrandedRegion)r;
			
			String seq = seqgen.execute(r).toUpperCase();
			if(seq.contains("N"))
				continue;
			for(int k=kmin; k<=kmax; k++){
				int currKStart = 0;
				for(int k_sub=kmin; k_sub<k; k_sub++){
					currKStart = currKStart  +  ((int)Math.pow(4, k_sub))*3;
				}
				// Count left kmers
				for(int i=0; i<((int)(seq.length()/3)-k+1); i++){
					String currK = seq.substring(i, i+k);
					String revCurrK =SequenceUtils.reverseComplement(currK);
					int  currKInt = RegionFileUtilities.seq2int(currK);
					int  revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					if(sr.getStrand() == '+')
						kmerCounts[currKStart+kmer]++;
					else
						kmerCounts[currKStart+ ((int)Math.pow(4, k))*2+kmer]++;
					
				}
				//Count central kmers
				for(int i=(int)(seq.length()/3); i<((int)(seq.length()*2/3)-k+1); i++){
					String currK = seq.substring(i, i+k);
					String revCurrK =SequenceUtils.reverseComplement(currK);
					int  currKInt = RegionFileUtilities.seq2int(currK);
					int  revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					kmerCounts[currKStart+ ((int)Math.pow(4, k))+kmer]++;
				}
				//Count right kmers
				for(int i=(int)(seq.length()*2/3); i<(seq.length()-k+1); i++){
					String currK = seq.substring(i, i+k);
					String revCurrK =SequenceUtils.reverseComplement(currK);
					int  currKInt = RegionFileUtilities.seq2int(currK);
					int  revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					if(sr.getStrand() == '+')
						kmerCounts[currKStart+ ((int)Math.pow(4, k))*2+kmer]++;
					else
						kmerCounts[currKStart+kmer]++;
				}
			}
			System.out.print(r.getLocationString());
			for(int i=0; i<kmerCounts.length; i++)
				System.out.print("\t"+kmerCounts[i]);
			System.out.println("");
			
			
			
		}
	}
		
	
	
	///////////////////////////////////////////////////////////////////////
	
	public void setWin(int w){window=w;}
	public void setReadLen(int rl){readLen=rl;}
	public void setFixedPB(int fpb){fixedPBCutoff = fpb;}
	public void setNormCounts(boolean x){normCounts = x;}

	//load positive
	public void loadPeaks(String fname, String motifHitsFile, boolean macs){
		if(fname == null && motifHitsFile != null){
			posSet=new ArrayList<Region>();
			posPeaks=new ArrayList<Point>();
			List<StrandedRegion> regs = RegionFileUtilities.loadStrandedRegionsFromMotifFile(gen, motifHitsFile, window);
			List<StrandedPoint> points = RegionFileUtilities.loadStrandedPointsFromMotifFile(gen, motifHitsFile, window);
			for(StrandedRegion r : regs)
				posSet.add(r);
			for(StrandedPoint p : points)
				posPeaks.add(p);
			posLines = RegionFileUtilities.loadLinesFromFile(motifHitsFile);
		}else{
			if(macs){
				posSet=new ArrayList<Region>();
				List<MACSPeakRegion> mprs = MACSParser.parseMACSOutput(fname,gen);
				for(MACSPeakRegion m : mprs){
					posSet.add(m.getPeak().expand(window/2));
				}
			}else{
				posSet = RegionFileUtilities.loadRegionsFromPeakFile(gen, fname, window);			
				posPeaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, fname, window);
			}
			posLines = RegionFileUtilities.loadLinesFromFile(fname);
		}	
	}//load negative
	public void loadNegs(String fname, boolean macs){
		if(macs){
			negSet=new ArrayList<Region>();
			List<MACSPeakRegion> mprs = MACSParser.parseMACSOutput(fname,gen);
			for(MACSPeakRegion m : mprs){
				negSet.add(m.getPeak().expand(window/2));
			}
		}else{
			negSet = RegionFileUtilities.loadRegionsFromPeakFile(gen, fname, window);			
			negPeaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, fname, window);
		}negLines = RegionFileUtilities.loadLinesFromFile(fname);	
	}//load towers
	public void loadTowers(String fname){
		towers = RegionFileUtilities.loadRegionsFromPeakFile(gen, fname, -1);
	}
	
	//load metapeak
	public double [] loadMetaPeakFile(String fname){
		double [] meta = new double [(window/metaBinSize)+1];
		try{
			File pFile = new File(fname);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        int lineCount=0, bCount=0;;
	        while ((line = reader.readLine()) != null) {
	            if(lineCount>=2 && bCount<=(window/metaBinSize)){
		        	line = line.trim();
		            String[] words = line.split("\\s+");
		            meta[bCount]=Double.valueOf(words[4]);
		            bCount++;
	            }		            
	            lineCount++;
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(meta);
	}
	
	
	//Print the read counts under each peak
	public void printXY(){
		try {
			FileWriter fout = new FileWriter("XYcounts.txt");
			for(Region r : posSet){
				//Initialize tower filter
				ArrayList<Region> currTowers = new ArrayList<Region>();
				for(Region t : towers)
					if(r.overlaps(t))
						currTowers.add(t);
				
				int IPcount=0, IP2count=0;
				for(ReadHit h : IP.loadHits(r)){
					boolean inTower=false;
					for(Region t : currTowers)
						if(h.overlaps(t))
							inTower=true;
					if(!inTower){IPcount++;}
				}
				for(ReadHit h : IP2.loadHits(r)){
					boolean inTower=false;
					for(Region t : currTowers)
						if(h.overlaps(t))
							inTower=true;
					if(!inTower){IP2count++;}
				}
				
				fout.write(String.format("%s\t%d\t%d\n", r, IPcount, IP2count));
			}fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	//Print the distribution of reads around the peaks
	//With tower & needle filters
	public void weightedCounts(double[] meta, String label, String outName){
		int topX=-1;
		boolean absoluteDist=true;
		boolean needleFiltering=true;
		int min=0, max=window/2;
		int numBins = ((max-min)/metaBinSize);
		int perBaseMax = getPoissonThreshold(Math.pow(10, -9), IP.getWeightTotal(), 1, gen.getGenomeLength(), 0.8, 1, 1);
		
		double [] positive = new double [numBins+1];
		for(int w=0; w<=numBins; w++){
			positive[w]=0;
		}
		try{
			FileWriter fw = new FileWriter(outName);
			fw.write("Region\t"+label+"\n");
			//Positive set
			int count =0;
			for(Region r : posSet){
				double score=0;
				if(topX==-1 || count<topX){
					Point p = posPeaks.get(count);
					List<ReadHit> hits = IP.loadHits(r);
					//Towers & Needles init
					ArrayList<Region> currTowers = new ArrayList<Region>();
					for(Region t : towers)
						if(r.overlaps(t))
							currTowers.add(t);
					int [] counts = new int[r.getWidth()+1];
					for(int i=0; i<=r.getWidth(); i++){counts[i]=0;}
					
					for(ReadHit h : hits){
						boolean inTower=false;
						for(Region t : currTowers)
							if(r.overlaps(t))
								inTower=true;
						if(!inTower){
							int offset=h.getFivePrime()-r.getStart()<0 ? 0 : h.getFivePrime()-r.getStart();
							if(offset>r.getWidth())
								offset=r.getWidth();
							counts[offset]++;
							if(!needleFiltering || (counts[offset] <= perBaseMax)){
								int dist = h.getStrand()=='+' ? p.getLocation()-h.getFivePrime() : h.getFivePrime()-p.getLocation();
								if(absoluteDist){dist = Math.abs(dist);}
								if(dist>=min && dist<=max){
									score+=meta[(dist-min)/metaBinSize];
								}
							}
						}
					}	
					fw.write(String.format("%s\t%.2f\n", p.getLocationString(),score));
				}count++;
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the distribution of reads around the peaks
	//With tower & needle filters
	public void counts(String label, String outName){
		int topX=-1;
		boolean needleFiltering=true;
		double ipWeightTotal=IP.getWeightTotal();
		int perBaseMax = fixedPBCutoff>0 ? fixedPBCutoff : getPoissonThreshold(Math.pow(10, -9), ipWeightTotal, 1, gen.getGenomeLength(), 0.8, 1, 1);
		try{
			FileWriter fw = new FileWriter(outName);
			fw.write("Region\t"+label+"\n");
			//Positive set
			int count =0;
			for(Region r : posSet){
				double score=0;
				if(topX==-1 || count<topX){
					Point p = posPeaks.get(count);
					List<ReadHit> hits = IP.loadHits(r);
					//Towers & Needles init
					ArrayList<Region> currTowers = new ArrayList<Region>();
					for(Region t : towers)
						if(r.overlaps(t))
							currTowers.add(t);
					int [] counts = new int[r.getWidth()+1];
					for(int i=0; i<=r.getWidth(); i++){counts[i]=0;}
					
					for(ReadHit h : hits){
						boolean inTower=false;
						for(Region t : currTowers)
							if(r.overlaps(t))
								inTower=true;
						if(!inTower){
							int offset=h.getFivePrime()-r.getStart()<0 ? 0 : h.getFivePrime()-r.getStart();
							if(offset>r.getWidth())
								offset=r.getWidth();
							counts[offset]++;
							if(!needleFiltering || (counts[offset] <= perBaseMax)){
								score++;
							}
						}
					}	
					if(normCounts)
						fw.write(String.format("%s\t%e\n", r.getLocationString(), score/ipWeightTotal));
					else
						fw.write(String.format("%s\t%.2f\n", r.getLocationString() ,score));
				}count++;
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Order the input peaks by read counts in surrounding window
	//With tower & needle filters
	public void orderByCounts(String label, String outName){
		int topX=-1;
		boolean needleFiltering=true;
		double ipWeightTotal=IP.getWeightTotal();
		int perBaseMax = fixedPBCutoff>0 ? fixedPBCutoff : getPoissonThreshold(Math.pow(10, -9), ipWeightTotal, 1, gen.getGenomeLength(), 0.8, 1, 1);
		ArrayList<Pair<Integer, Double>> peakScores = new ArrayList<Pair<Integer, Double>>();
		try{
			FileWriter fw = new FileWriter(outName);
			//Positive set
			int count =0;
			for(Region r : posSet){
				double score=0;
				if(topX==-1 || count<topX){
					Point p = posPeaks.get(count);
					List<ReadHit> hits = IP.loadHits(r);
					//Towers & Needles init
					ArrayList<Region> currTowers = new ArrayList<Region>();
					for(Region t : towers)
						if(r.overlaps(t))
							currTowers.add(t);
					int [] counts = new int[r.getWidth()+1];
					for(int i=0; i<=r.getWidth(); i++){counts[i]=0;}
					
					for(ReadHit h : hits){
						boolean inTower=false;
						for(Region t : currTowers)
							if(r.overlaps(t))
								inTower=true;
						if(!inTower){
							int offset=h.getFivePrime()-r.getStart()<0 ? 0 : h.getFivePrime()-r.getStart();
							if(offset>r.getWidth())
								offset=r.getWidth();
							counts[offset]++;
							if(!needleFiltering || (counts[offset] <= perBaseMax)){
								score++;
							}
						}
					}
					peakScores.add(new Pair<Integer,Double>(count, score));
				}count++;
			}
			//Sort
			Collections.sort(peakScores, new Comparator<Pair<Integer,Double>>(){
				public int compare(Pair<Integer,Double> o1, Pair<Integer,Double> o2){return o2.cdr().compareTo(o1.cdr());}
			});
				
			//Print
			for(int x=0; x<peakScores.size(); x++){
				fw.write(posLines.get(peakScores.get(x).car())+"\n");
			}
			
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the distribution of reads around the peaks
	//NO tower & needle filters -- built for speed
	public void quickcounts(String label, String outName){
		int topX=-1;
		double ipWeightTotal=IP.getWeightTotal();
		ArrayList<Pair<String, Double>> results = new ArrayList<Pair<String, Double>>(); 
		try{
			
			//Positive set
			int count =0;
			for(Region r : posSet){
				double score=0;
				if(topX==-1 || count<topX){
					Point p = posPeaks.get(count);
					score = IP.sumWeights(r);	
					if(normCounts)
						results.add(new Pair<String,Double>(r.getLocationString(),score/ipWeightTotal));
					else
						results.add(new Pair<String,Double>(r.getLocationString(),score));
				}count++;
			}
			FileWriter fw = new FileWriter(outName);
			fw.write("Region\t"+label+"\n");
			for(Pair<String,Double> res:results){
				if(normCounts)
					fw.write(String.format("%s\t%e\n", res.car(), res.cdr()));
				else
					fw.write(String.format("%s\t%.2f\n", res.car(), res.cdr()));
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the distribution of reads around the peaks
	//With tower & needle filters
	public void buildMetaPeak(String outName){
		System.out.println("Building Meta-peak");
		int topX=-1;
		boolean absoluteDist=true;
		boolean needleFiltering=true;
		int min=0, max=window/2;
		int numBins = ((max-min)/metaBinSize);
		int perBaseMax = getPoissonThreshold(Math.pow(10, -9), IP.getWeightTotal(), 1, gen.getGenomeLength(), 0.8, 1, 1);
		
		//Initialize: add pseudocount of 1 to everything
		double [] positive = new double [numBins+1];
		double [] negative = new double [numBins+1];
		for(int w=0; w<=numBins; w++){
			positive[w]=1; negative[w]=1;
		}double posTotalHits=numBins+1, posTotalPeaks=0, negTotalHits=numBins+1,  negTotalPeaks=0;
		
		//Positive set
		int count =0;
		for(Region r : posSet){
		    //			System.out.println("POSITIVE: "+r);
			if(topX==-1 || count<topX){
				Point p = posPeaks.get(count);
				List<ReadHit> hits = IP.loadHits(r);
				//Towers & Needles init
				ArrayList<Region> currTowers = new ArrayList<Region>();
				for(Region t : towers)
					if(r.overlaps(t))
						currTowers.add(t);
				int [] counts = new int[r.getWidth()+1];
				for(int i=0; i<=r.getWidth(); i++){counts[i]=0;}
				
				for(ReadHit h : hits){
					boolean inTower=false;
					for(Region t : currTowers)
						if(r.overlaps(t))
							inTower=true;
					if(!inTower){
						int offset=h.getFivePrime()-r.getStart()<0 ? 0 : h.getFivePrime()-r.getStart();
						if(offset>r.getWidth())
							offset=r.getWidth();
						counts[offset]++;
						if(!needleFiltering || (counts[offset] <= perBaseMax)){
							int dist = h.getStrand()=='+' ? p.getLocation()-h.getFivePrime() : h.getFivePrime()-p.getLocation();
							if(absoluteDist){dist = Math.abs(dist);}
							if(dist>=min && dist<=max){
								positive[(dist-min)/metaBinSize]++;
								posTotalHits++;
							}
						}
					}
				}posTotalPeaks++;			
			}count++;
		}
		//Negative set
		count =0;
		for(Region r : negSet){
		    //System.out.println("NEGATIVE: "+r);
			if(topX==-1 || count<topX){
				Point p = negPeaks.get(count);
				List<ReadHit> hits = IP.loadHits(r);
				//Towers & Needles init
				ArrayList<Region> currTowers = new ArrayList<Region>();
				for(Region t : towers)
					if(r.overlaps(t))
						currTowers.add(t);
				int [] counts = new int[r.getWidth()+1];
				for(int i=0; i<=r.getWidth(); i++){counts[i]=0;}
				
				for(ReadHit h : hits){
					boolean inTower=false;
					for(Region t : currTowers)
						if(r.overlaps(t))
							inTower=true;
					if(!inTower){
						int offset=h.getFivePrime()-r.getStart()<0 ? 0 : h.getFivePrime()-r.getStart();
						if(offset>r.getWidth())
							offset=r.getWidth();
						counts[offset]++;
						if(!needleFiltering || (counts[offset] <= perBaseMax)){
							int dist = h.getStrand()=='+' ? p.getLocation()-h.getFivePrime() : h.getFivePrime()-p.getLocation();
							if(absoluteDist){dist = Math.abs(dist);}
							if(dist>=min && dist<=max){
								negative[(dist-min)/metaBinSize]++;
								negTotalHits++;
							}
						}
					}
				}negTotalPeaks++;
			}count++;
		}
		
		//Normalize
		for(int w=0; w<=numBins; w++){
			positive[w]=positive[w]/posTotalPeaks; 
			negative[w]=negative[w]/negTotalPeaks;
		}
		
		try{
			FileWriter fw = new FileWriter(outName);
			fw.write("Tag distributions around peaks\n");
			fw.write("RelPos\tPositive\tNegative\tDifference\tOverRep\n");
			for(int w=0; w<numBins; w++){
				int rel = (w*metaBinSize)+min;
				fw.write(rel+"\t"+positive[w]+"\t"+negative[w]+"\t"+(positive[w]-negative[w])+"\t"+(positive[w]/negative[w])+"\n");
			}
			fw.close();
			System.out.println("Results written to "+outName);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected void loadExperiments(String [] args){
		List<SeqLocator> dbexpts = Args.parseSeqExpt(args,"dbexpt");
        List<SeqLocator> dbexpts2 = Args.parseSeqExpt(args,"dbexpt2");
        List<SeqLocator> rdbexpts = Args.parseSeqExpt(args,"rdbexpt");
        List<SeqLocator> rdbexpts2 = Args.parseSeqExpt(args,"rdbexpt2");
        List<File> expts = Args.parseFileHandles(args, "expt");
        List<File> expts2 = Args.parseFileHandles(args, "expt2");
        String fileFormat = Args.parseString(args, "format", "ELAND");
        if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
        	IP = new DeepSeqExpt(gen, expts, false, fileFormat, readLen);
        }else if(dbexpts.size()>0 && expts.size() == 0){
        	IP = new DeepSeqExpt(gen, dbexpts, "db", readLen);
        	dbconnected=true;
        }else if(rdbexpts.size()>0 && expts.size() == 0){
        	IP = new DeepSeqExpt(gen, rdbexpts, "readdb", readLen);
        	dbconnected=true;
        }else{}
        if(expts2.size()>0 && dbexpts2.size() == 0){
        	IP2 = new DeepSeqExpt(gen, expts2, false, fileFormat, readLen); twoExpts=true;
        }else if(dbexpts2.size()>0 && expts2.size() == 0){
        	IP2 = new DeepSeqExpt(gen, dbexpts2, "db", readLen); twoExpts=true;
        	dbconnected=true;
        }else if(rdbexpts2.size()>0 && expts2.size() == 0){
        	IP2 = new DeepSeqExpt(gen, rdbexpts2, "readdb", readLen); twoExpts=true;
        	dbconnected=true;
        }else{
        	if(dbexpts2.size()>0 && expts2.size()>0){
        		System.err.println("Cannot mix files and db loading yet...");System.exit(1);
        	}else{
        		twoExpts=false; IP2=null;
        	}
        }
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
	
	//Clean up the loaders
	public void close(){
		if(IP!=null)
			IP.closeLoaders();
		if(IP2!=null)
			IP2.closeLoaders();
	}
}
