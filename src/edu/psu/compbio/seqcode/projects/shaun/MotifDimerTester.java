package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.FASTALoader;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.motifs.BackgroundModelIO;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class MotifDimerTester {
	private Organism org;
	private Genome gen;
	private WeightMatrix motifA;
	private WeightMatrix motifB;
	private ArrayList<Region> regions;
	private ArrayList<String> seqs = new ArrayList<String>();
        private ArrayList<Region> negregions = new ArrayList<Region>();
	private ArrayList<String> negseqs = new ArrayList<String>();
	private ArrayList<String> inFileLines;
	private double motifThresA=-100000;
	private double motifThresB=-100000;
	private int averageRegWidth=200, maxRegWidth=0; 
	private int maxSpacer = 9;
	private int numRandom = 100000;
	private SequenceGenerator seqgen = new SequenceGenerator();
	private int [] DRcount, ERcount, IRcount;
	private int [] DRpeaks, ERpeaks, IRpeaks;
	private boolean printMatches=false;
	private boolean printSequences=false;
	private boolean havePeaks = false;
	private boolean makeRandom=false;
	private boolean noOverlappingDimers = false;
	private int win =200;
	private ArrayList<DimerHit> totalDimers = new ArrayList<DimerHit>();
	private ArrayList<DimerHit> finalDimers = new ArrayList<DimerHit>();
	
	public static void main(String[] args) throws IOException, ParseException {
		MotifDimerTester tester;
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||(!ap.hasKey("motifnameA")&&!ap.hasKey("motiffileA")) || (!ap.hasKey("motifnameB")&&!ap.hasKey("motiffileB")) || (!ap.hasKey("peaks")&&!ap.hasKey("scoreseqs")&&!ap.hasKey("random"))) { 
            System.err.println("Usage:\n " +
                               "MotifThresholdFinder " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--motiffileA <file with motifs> AND --back <background>\nOR\n"+
                               "--motifnameA <weightmatrix name> "+
                               "--motifversionA <weightmatrix version> \n"+
                               "--motiffileB <file with motifs> AND --back <background>\nOR\n"+
                               "--motifnameB <weightmatrix name> "+
                               "--motifversionB <weightmatrix version> \n"+
                               "--peaks <file containing coordinates of peakst> \n"+
                               "--win <window of sequence to take around peaks> \n" +
                               "--scoreseqs <fasta file to score> \n" +
                               "--random <number of random positions> \n" +
                               "--printmatches [flag to print motif match IDs] \n"+
                               "--printsequences [flag to print motif match sequences] \n" +
                               "--nooverlapping [flag to discount overlapping dimers] \n"+
                               "--minscoreA <minimum motif score> \n" +
                               "--minscoreB <minimum motif score> \n");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String motifnameA=null, motifnameB=null;
        String motifversionA=null;
        String motifversionB=null;
        String motiffileA=null; String motiffileB=null;
        String backfile=null;
        boolean loadFromFileA=false, loadFromFileB=false;
        String posFile = (ap.hasKey("peaks")) ? ap.getKeyValue("peaks") : null;
        boolean makeRandom= ap.hasKey("random");
        int numRand = ap.hasKey("random") ? new Integer(ap.getKeyValue("random")).intValue():-1;
        boolean printSequences= ap.hasKey("printsequences");
        boolean printMatches= ap.hasKey("printmatches");
        boolean usingWin= ap.hasKey("win");
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        double minscoreA = (ap.hasKey("minscoreA")) ? new Double(ap.getKeyValue("minscoreA")).doubleValue() : -10000;
        double minscoreB = (ap.hasKey("minscoreB")) ? new Double(ap.getKeyValue("minscoreB")).doubleValue() : -10000;
        boolean scoreSeqs = ap.hasKey("scoreseqs");
        boolean overlaps = !ap.hasKey("nooverlapping");
        String seqFile=null;
        if(scoreSeqs)
        	seqFile = ap.getKeyValue("scoreseqs");
        if(ap.hasKey("motifnameA")){
	        if(ap.hasKey("motifversionA")){motifversionA = ap.getKeyValue("motifversionA");}
	        motifnameA = ap.getKeyValue("motifnameA");
	        if (motifnameA.indexOf(';') != -1) {
	            String[] pieces = motifnameA.split(";");
	            motifnameA = pieces[0];
	            motifversionA = pieces[1];
	        }
        }else if(ap.hasKey("motiffileA")){
	       	loadFromFileA=true;
	       	motiffileA = ap.getKeyValue("motiffileA");
	       	backfile = ap.getKeyValue("back");
	    }
        if(ap.hasKey("motifnameB")){
	        if(ap.hasKey("motifversionB")){motifversionB = ap.getKeyValue("motifversionB");}
	        motifnameB = ap.getKeyValue("motifnameB");
	        if (motifnameB.indexOf(';') != -1) {
	            String[] pieces = motifnameB.split(";");
	            motifnameB = pieces[0];
	            motifversionB = pieces[1];
	        }
        }else if(ap.hasKey("motiffileB")){
	       	loadFromFileB=true;
	       	motiffileB = ap.getKeyValue("motiffileB");
	       	backfile = ap.getKeyValue("back");
	    }
        
        
        ArrayList<Region> posRegs = new ArrayList<Region>();
        ArrayList<String> lines= new ArrayList<String>();;
		try {
			//Load genome
			Organism currorg = Organism.getOrganism(species);
			Genome currgen = currorg.getGenome(genome);
			
			//Load motifs
			WeightMatrix matrixA = null;
			WeightMatrix matrixB = null;
			if(loadFromFileA){
				matrixA = loadMotifFromFile(motiffileA, backfile, currgen);
			}else{			
				int wmidA = WeightMatrix.getWeightMatrixID(currorg.getDBID(), motifnameA, motifversionA);
		        matrixA = WeightMatrix.getWeightMatrix(wmidA);
			}
	        if(loadFromFileB){
	        	matrixB = loadMotifFromFile(motiffileB, backfile, currgen);
	        }else{
	        	int wmidB = WeightMatrix.getWeightMatrixID(currorg.getDBID(), motifnameB, motifversionB);
	        	matrixB = WeightMatrix.getWeightMatrix(wmidB);
	        }
        
	        //Initialize
	        tester = new MotifDimerTester(currorg, currgen, matrixA, matrixB, minscoreA, minscoreB);
	        tester.setMakeRandom(makeRandom);
	        tester.setPrintMatches(printMatches);
	        tester.setPrintSequences(printSequences);
	        tester.setNoOverlaps(!overlaps);
	        if(usingWin)
	        	tester.setWin(win);
	        if(makeRandom)
	            tester.setNumRandom(numRand);
	        
	        if(posFile!=null){
	            //Load the positive hits
		        File pFile = new File(posFile);
				if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
		        BufferedReader reader = new BufferedReader(new FileReader(pFile));
		        String line= reader.readLine();
		        while ((line = reader.readLine()) != null) {
		            line = line.trim();
		            lines.add(line);
		            String[] words = line.split("\\s+");
		            if(usingWin){
		            	PointParser pparser = new PointParser(currgen);
		            	Point p = pparser.execute(words[2]);
		                int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
		            	int rend = p.getLocation()+(win/2)>currgen.getChromLength(p.getChrom()) ? currgen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
		            	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
		            	posRegs.add(r);
		            }else{
		            	RegionParser parser = new RegionParser(currgen);
		            	Region r = parser.execute(words[0]);
		            	if(r!=null){posRegs.add(r);}
		            	Point p = r.getMidpoint();
		            }
		        }
		        tester.loadPeaks(lines, posRegs);
		        tester.execute();
	        }else if(scoreSeqs && seqFile!=null){
	        	ArrayList<String> s = new ArrayList<String>();
	        	FASTALoader loader = new FASTALoader();
	    		File f = new File(seqFile);
	    		Iterator<Pair<String, String>> it = loader.execute(f);
	    		while(it.hasNext()){
	    			Pair<String, String> p = it.next();
	    			String name = p.car();
	    			String seq = p.cdr();
	    			s.add(seq);
	    		}
	        	tester.execute(s);	
	        }else if(makeRandom){
	            tester.execute();
	        }   
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	//Constructor 
	public MotifDimerTester(Organism o, Genome g, WeightMatrix wmA, WeightMatrix wmB, double minA, double minB){
		org=o;
		gen=g;
		motifA= wmA;
		motifB= wmB;
		motifThresA = minA;
		motifThresB = minB;
	}
	
	public void loadPeaks(ArrayList<String> inL,  ArrayList<Region> posreg){
		inFileLines = inL;
		regions=posreg;
		averageRegWidth = findAvgWidth(regions);
		maxRegWidth = findMaxWidth(regions);
		havePeaks=true;
	}
	public void setHavePeaks(boolean p){havePeaks=p;}
	public void setPrintMatches(boolean p){printMatches=p;}
	public void setPrintSequences(boolean p){printSequences=p;}
	public void setNoOverlaps(boolean no){noOverlappingDimers = no;}
	public void setMakeRandom(boolean p){makeRandom=p;}
	public void setNumRandom(int nr){numRandom=nr;}
	public void setWin(int w){win=w;}
	
	private int findAvgWidth(ArrayList<Region> regs){
		long total=0, num=0;
		for(Region r: regs){
			total+=r.getWidth();
			num++;
		}
		return((int)(total/num));
	}
	private int findMaxWidth(ArrayList<Region> regs){
		int max=0;
		for(Region r: regs)
			if(r.getWidth()>max){max = r.getWidth();}
		return(max);
	}
	
	//execute a dimer search
	public void execute(){
		String oString = " (overlapping dimers allowed)";
		if(noOverlappingDimers)
			oString =  " (no overlapping dimers allowed)";
		
		if(havePeaks){
			System.out.println("Dimers in "+regions.size()+" user-supplied regions"+oString);
			seqs = new ArrayList<String>();
			for(int i=0; i<regions.size(); i++){
				Region r = regions.get(i);
				seqs.add(seqgen.execute(r));
			}
			countDimerHits(seqs);
		}
		if(makeRandom){
			//Make random regions & test
			System.out.println("\n\nDimers in "+numRandom+" Random Regions"+oString);
			negregions= randomRegionPick(numRandom, averageRegWidth);
			negseqs = new ArrayList<String>();
			for(int i=0; i<negregions.size(); i++){
				Region r = negregions.get(i);
				negseqs.add(seqgen.execute(r));
			}
			countDimerHits(negseqs);
		}
	}
	
	//execute a dimer search
	public void execute(ArrayList<String> s){
		String oString = " (overlapping dimers allowed)";
		if(noOverlappingDimers)
			oString =  " (no overlapping dimers allowed)";
		
		System.out.println("Dimers in "+s.size()+" user-supplied sequences"+oString);
		seqs = s;
		countDimerHits(s);
		if(makeRandom){
			//Make random regions & test
			System.out.println("\n\nDimers in "+numRandom+" Random Regions"+oString);
			negregions= randomRegionPick(numRandom, averageRegWidth);
			negseqs = new ArrayList<String>();
			for(int i=0; i<negregions.size(); i++){
				Region r = negregions.get(i);
				negseqs.add(seqgen.execute(r));
			}
			countDimerHits(negseqs);
		}
	}
	//Simple count of hits
	public void countDimerHits(ArrayList<String> sequences){
		int monoCountA=0, monoCountB=0;
		int peaksWithA=0, peaksWithB=0;
		DRcount = new int [maxSpacer+1]; ERcount = new int [maxSpacer+1]; IRcount = new int [maxSpacer+1];
		DRpeaks = new int [maxSpacer+1]; ERpeaks = new int [maxSpacer+1]; IRpeaks = new int [maxSpacer+1];
		for(int s=0; s<=maxSpacer; s++){
			DRcount[s]=0; ERcount[s]=0; IRcount[s]=0;
			DRpeaks[s]=0; ERpeaks[s]=0; IRpeaks[s]=0;
		}
		
		int dimerHits=0, peaksWithDimerHits=0, totalPeaks=0;
		WeightMatrixScorer scorerA = new WeightMatrixScorer(motifA);
		WeightMatrixScorer scorerB = new WeightMatrixScorer(motifB);
		
		for(int i=0; i<sequences.size(); i++){
			totalPeaks++;
			String seq = sequences.get(i);
			double [] motifAHits = new double [seq.length()+1];
			double [] motifBHits = new double [seq.length()+1];
			for(int z=0; z<seq.length(); z++){
				motifAHits[z]=0; motifBHits[z]=0;
			}
			
			WeightMatrixScoreProfile profilerA = scorerA.execute(seq);
			boolean monoA=false;
			for(int z=0; z<seq.length(); z++){
				double currScore= profilerA.getMaxScore(z);
				if(currScore>motifThresA){
					if(profilerA.getMaxStrand(z)=='+')
						motifAHits[z]=1;
					else
						motifAHits[z]=-1;
					monoCountA++;
					monoA=true;
				}
			}if(monoA){peaksWithA++;}
			
			WeightMatrixScoreProfile profilerB = scorerB.execute(seq);
			boolean monoB=false;
			for(int z=0; z<seq.length(); z++){
				double currScore= profilerB.getMaxScore(z);
				if(currScore>motifThresB){
					if(profilerB.getMaxStrand(z)=='+')
						motifBHits[z]=1;
					else
						motifBHits[z]=-1;
					monoCountB++;
					monoB=true;
				}
			}if(monoB){peaksWithB++;}	
			
			ArrayList<String> thisSequenceDimers=new ArrayList<String>(); //To avoid duplicates, this container will hold IDs based on left-most coord
			int x,y; boolean goodMotif=false;
			//Check for dimers
			for(x=0; x<seq.length(); x++){				
				if(motifAHits[x]==1 || motifAHits[x]==-1){
					//Direct-repeat dimers
					// F1 :: F2  &  R1 :: R2
				        // -->  -->     <--  <--
					for(y=x+motifA.length(); y<=x+motifA.length()+maxSpacer; y++){if(y>=0 && y<seq.length()){
						if(motifBHits[y]==motifAHits[x]){
							int gap = y-(x+motifA.length());
							String ID = "DR"+gap+"_i"+x;
							if(!thisSequenceDimers.contains(ID)){
								thisSequenceDimers.add(ID);
								goodMotif=true; 
								String subseq = seq.substring(x, y+motifB.length());
								char strand = '+';
								if(motifAHits[x]!=1){
									String revseq = SequenceUtils.reverseComplement(subseq);
									subseq = revseq;
									strand='-';
								}
								totalDimers.add(new DimerHit(i, x, y+motifB.length(),strand, profilerA.getMaxScore(x)+profilerA.getMaxScore(y), 'D', gap, subseq));
							}
						}
					}}
					// F2 :: F1  &  R2 :: R1
					// -->  -->     <--  <--
					for(y=x-(motifB.length()+maxSpacer); y<=x-motifB.length(); y++){if(y>=0 && y<seq.length()){
						if(motifBHits[y]==motifAHits[x]){
							int gap =-1*(y-(x-motifB.length())); 
							String ID = "DR"+gap+"_i"+y;
							if(!thisSequenceDimers.contains(ID)){
								thisSequenceDimers.add(ID);
								goodMotif=true; 
								String subseq = seq.substring(y, x+motifA.length());
								char strand = '+';
								if(motifAHits[x]!=1){
									String revseq = SequenceUtils.reverseComplement(subseq);
									subseq = revseq;
									strand='-';
								}
								totalDimers.add(new DimerHit(i, y, x+motifA.length(),strand, profilerA.getMaxScore(x)+profilerA.getMaxScore(y), 'D', gap, subseq));
							}
						}						
					}}
					
					//Repeat inverted (incl palindromes)
					// F1 :: R2  &  R1 :: F2
					// -->  <--     <--  -->
					for(y=x+motifA.length(); y<=x+motifA.length()+maxSpacer; y++){if(y>=0 && y<seq.length()){
						if(motifAHits[x]==1 && motifBHits[y]==-1){
							int gap = y-(x+motifA.length());
							String ID = "ER"+gap+"_i"+x;
							if(!thisSequenceDimers.contains(ID)){
								thisSequenceDimers.add(ID);
								goodMotif=true; 
								String subseq = seq.substring(x, y+motifB.length());
								char strand = '+';
								if(motifAHits[x]!=1){
									String revseq = SequenceUtils.reverseComplement(subseq);
									subseq = revseq;
									strand='-';
								}
								totalDimers.add(new DimerHit(i, x, y+motifB.length(),strand, profilerA.getMaxScore(x)+profilerA.getMaxScore(y), 'E', gap, subseq));
							}
						}else if(motifAHits[x]==-1 && motifBHits[y]==1){
							int gap = y-(x+motifA.length());
							String ID = "IR"+gap+"_i"+x;
							if(!thisSequenceDimers.contains(ID)){
								thisSequenceDimers.add(ID);
								goodMotif=true; 
								String subseq = seq.substring(x, y+motifB.length());
								char strand = '+';
								if(motifAHits[x]!=1){
									String revseq = SequenceUtils.reverseComplement(subseq);
									subseq = revseq;
									strand='-';
								}
								totalDimers.add(new DimerHit(i, x, y+motifB.length(),strand, profilerA.getMaxScore(x)+profilerA.getMaxScore(y), 'I', gap, subseq));
							}
						}
					}}
					// F2 :: R1  &  R2 :: F1
					// -->  <--     <--  -->
					for(y=x-(motifB.length()+maxSpacer); y<=x-motifB.length(); y++){if(y>=0 && y<seq.length()){
						if(motifAHits[x]==-1 && motifBHits[y]==1){
							int gap = -1*(y-(x-motifB.length()));
							String ID = "ER"+gap+"_i"+y;
							if(!thisSequenceDimers.contains(ID)){
								thisSequenceDimers.add(ID);
								goodMotif=true; 
								String subseq = seq.substring(y, x+motifA.length());
								char strand = '+';
								if(motifAHits[x]!=1){
									String revseq = SequenceUtils.reverseComplement(subseq);
									subseq = revseq;
									strand='-';
								}
								totalDimers.add(new DimerHit(i, y, x+motifA.length(),strand, profilerA.getMaxScore(x)+profilerA.getMaxScore(y), 'E', gap, subseq));
							}
						}else if(motifAHits[x]==1 && motifBHits[y]==-1){
							int gap = -1*(y-(x-motifB.length()));
							String ID = "IR"+gap+"_i"+y;
							if(!thisSequenceDimers.contains(ID)){
								thisSequenceDimers.add(ID);
								goodMotif=true; 
								String subseq = seq.substring(y, x+motifA.length());
								char strand = '+';
								if(motifAHits[x]!=1){
									String revseq = SequenceUtils.reverseComplement(subseq);
									subseq = revseq;
									strand='-';
								}
								totalDimers.add(new DimerHit(i, y, x+motifA.length(),strand, profilerA.getMaxScore(x)+profilerA.getMaxScore(y), 'I', gap, subseq));
							}
						}
					}}
				}
			}
			if(goodMotif){peaksWithDimerHits++;}
		}//totalDimers populated
		
		//Filter if necessary
		if(noOverlappingDimers){
			//Process all dimers
			while(totalDimers.size()>0){
				//Find the max score
				double maxScore =-10000000; 
				int maxIndex=0;
				for(int d=0; d<totalDimers.size(); d++){
					if(totalDimers.get(d).score>maxScore){
						maxScore = totalDimers.get(d).score; maxIndex=d;
					}
				}
				DimerHit max = totalDimers.get(maxIndex);
				//Check if max overlaps any previously added dimer
				boolean overlaps=false;
				for(DimerHit e : finalDimers){
					if(max.overlaps(e)){
						overlaps=true;
						break;
					}
				}
				//Add/discard
				if(!overlaps)
					finalDimers.add(max);
				totalDimers.remove(maxIndex);
			}
		}else{
			finalDimers = totalDimers;
		}
		dimerHits = finalDimers.size();
		
		//Count the remaining matches and positive peaks
		for(DimerHit d : finalDimers){
			if(d.type == 'D')				
				DRcount[d.spacer]++;
			else if(d.type == 'I')
				IRcount[d.spacer]++;
			else if(d.type == 'E')
				ERcount[d.spacer]++;
		}
		char [] types = {'D', 'I', 'E'};
		for(char t : types){
			for(int s=0; s<=maxSpacer; s++){
				HashMap<Integer, Boolean> positives= new HashMap<Integer, Boolean>();
				for(DimerHit d : finalDimers){
					if(d.type==t && d.spacer==s){
						positives.put(d.seqID, true);
					}
				}
				if(t=='D')
					DRpeaks[s]= positives.keySet().size();
				if(t=='I')
					IRpeaks[s]= positives.keySet().size();
				if(t=='E')
					ERpeaks[s]= positives.keySet().size();
			}
		}
		
		
		System.out.println(monoCountA +" monomer hits in "+peaksWithA+" peaks for "+motifA.getName()+" with threshold = "+motifThresA);
		System.out.println(monoCountB +" monomer hits in "+peaksWithB+" peaks for "+motifB.getName()+" with threshold = "+motifThresB);
		System.out.println(dimerHits+" dimer hits in "+peaksWithDimerHits+" sequences from "+totalPeaks+" total sequences.");
		System.out.println("\nType\thits\thitsFreq\tpeaks\tpeaksFreq");
		for(int s=0; s<=maxSpacer; s++){
			double cavg = (double)DRcount[s]/(double)totalPeaks;
			double pavg = (double)DRpeaks[s]/(double)totalPeaks;
			System.out.println("DR"+s+"\t"+DRcount[s]+"\t"+cavg+"\t"+DRpeaks[s]+"\t"+pavg);
		}for(int s=0; s<=maxSpacer; s++){
			double cavg = (double)IRcount[s]/(double)totalPeaks;
			double pavg = (double)IRpeaks[s]/(double)totalPeaks;
			System.out.println("IR"+s+"\t"+IRcount[s]+"\t"+cavg+"\t"+IRpeaks[s]+"\t"+pavg);
		}for(int s=0; s<=maxSpacer; s++){
			double cavg = (double)ERcount[s]/(double)totalPeaks;
			double pavg = (double)ERpeaks[s]/(double)totalPeaks;
			System.out.println("ER"+s+"\t"+ERcount[s]+"\t"+cavg+"\t"+ERpeaks[s]+"\t"+pavg);
		}
		
		//Sequences
		if(printSequences){
			for(char t : types){
				for(int s=0; s<=maxSpacer; s++){
					int size=0;
					if(t=='D')
						size = DRcount[s];
					if(t=='I')
						size = IRcount[s];
					if(t=='E')
						size = ERcount[s];
					if(size>0){
						System.out.println("\n>"+t+"R"+s);
						for(DimerHit d : finalDimers){
							if(d.type==t && d.spacer==s){
								System.out.println(d.seq);
							}
						}
					}
				}
			}
		}
	}
	
	//Randomly pick a set of Regions
	private ArrayList<Region> randomRegionPick(int numSamples, int sampleSize){
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
						
						//is this region in the blacklist 
						boolean valid=true;
						if(potential.getEnd()>=chromoSize[c]){valid=false;}
						
						//or already picked?
						/*for(Region r : regs){
							if(potential.overlaps(r)){valid=false;}
						}*/
						
						if(valid){
							validSamples++;
							regs.add(potential);
							//System.out.println(potential.getChrom()+":"+potential.getStart()+"-"+potential.getEnd());
						}
					}
				}total+=chromoSize[c];
			}				
		}
		return(regs);
	}
	//Load freq matrix
	public static WeightMatrix loadMotifFromFile(String filename, String backname, Genome gen) throws IOException, ParseException {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		MarkovBackgroundModel back = BackgroundModelIO.parseMarkovBackgroundModel(backname, gen);
    	
    	motifImport.setBackground(back);
		return(motifImport.readTransfacMatrices(filename).getFirst());		
	}
	
	//A dimer hit
	public class DimerHit{
		public int seqID;
		public int start;
		public int end;
		public char strand;
		public double score;
		public char type; //D=direct, I=inverted, E=everted repeat
		public int spacer;
		public String seq;
		
		public DimerHit(int seqID, int start, int end, char strand, double score, char type, int spacer, String seq){
			this.seqID=seqID;
			this.start=start;
			this.end=end;
			this.strand=strand;
			this.score=score;
			this.type = type;
			this.spacer=spacer;
			this.seq = seq;
		}
		
		public boolean overlaps(DimerHit d) {
		    if (seqID != d.seqID)
		    	return false;
		    if (start <= d.start && end >= d.start)
		      return true;
		    if (d.start <= start && d.end >= start)
		      return true;
		    return false;
		}
	}
}
