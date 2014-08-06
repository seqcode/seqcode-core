package edu.psu.compbio.seqcode.projects.shaun;

import java.awt.Container;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.genome.Genome;
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
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;

public class MotifThresholdFinder {
	private WeightMatrix motif = null;
	private Genome gen=null;
	private ArrayList<String> seqSet = new ArrayList<String>(); 
	private int window=100;
    private int stepSize = 10;
	public ArrayList<Region> posSet = new ArrayList<Region>();
    	public ArrayList<Region> negSet = new ArrayList<Region>();
    public ArrayList<String> posSeq = new ArrayList<String>();
    public ArrayList<String> negSeq = new ArrayList<String>();
	private boolean ROC=false;
	
	public static void main(String[] args) throws IOException, ParseException {
		MotifThresholdFinder finder;
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||(!ap.hasKey("motifname")&&!ap.hasKey("motiffile"))) { 
            System.err.println("Usage:\n " +
                               "MarkovMotifThresholdFinder \n" +
                               "--species <organism name> \n" +
                               "--genome <genome version> \n"+
                               "--motifname <weightmatrix name> \n"+
                               "--motifversion <weightmatrix version> \n" +
                               "--motiffile <file containing motifs> \n"+
                               "--back <background Markov model> \n"+
                               "--win <window of sequence around positive/negative points> \n"+
                               "--positive <file containing coordinates of positive set> \n"+
                               "--negative <file containing coordinates of negative set (optional)> \n"+
                               "--numrand <number of random points to sample> \n" +
                               "--printroc \n");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String motifversion=null;
        if(ap.hasKey("motifversion")){motifversion = ap.getKeyValue("motifversion");}
        String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        int numSim = ap.hasKey("numrand") ? new Integer(ap.getKeyValue("numrand")).intValue() : 1000000;
        boolean printROC= ap.hasKey("printroc");
        boolean loadFromFile = ap.hasKey("motiffile");
        boolean posHitsGiven = ap.hasKey("positive");
        String posFile = ap.hasKey("positive") ? ap.getKeyValue("positive"):null;
        boolean negHitsGiven = ap.hasKey("negative");
        String negFile = ap.hasKey("negative") ? ap.getKeyValue("negative"):null;
        
        try {
			//Load genome
			Organism currorg = Organism.getOrganism(species);
			Genome currgen = currorg.getGenome(genome);

	        //Load the background model
	        MarkovBackgroundModel backMod;
	        if(backFile == null){
	          backMod = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(Organism.findGenome(genome)));
	        }else{
	        	backMod = BackgroundModelIO.parseMarkovBackgroundModel(backFile, Organism.findGenome(genome));
	        }
	        
	        //Load datasets
	        ArrayList<Region> pos = new ArrayList<Region>();
	    	ArrayList<Region> neg = new ArrayList<Region>();
	    	if(posHitsGiven)
	    		pos = loadRegionsFromPeakFile(currgen, posFile, win);
	        if(negHitsGiven)
	        	neg = loadRegionsFromPeakFile(currgen, negFile, win);
	        else{
	        	neg = randomRegionPick(currgen, numSim, win);
	        }
	        
			//Load motifs
	        List<WeightMatrix> motifList=new ArrayList<WeightMatrix>();
	        if(loadFromFile){
	        	String motifFile = ap.getKeyValue("motiffile");
	        	FreqMatrixImport motifImport = new FreqMatrixImport();
	        	motifImport.setBackground(backMod);
	    		motifList.addAll(motifImport.readTransfacMatrices(motifFile));
	    		
	        }else{
		        String motifname = ap.getKeyValue("motifname");
		        if (motifname.indexOf(';') != -1) {
		            String[] pieces = motifname.split(";");
		            motifname = pieces[0];
		            motifversion = pieces[1];
		        }
				int wmid = WeightMatrix.getWeightMatrixID(currorg.getDBID(), motifname, motifversion);
		        motifList.add(WeightMatrix.getWeightMatrix(wmid));
	        }
	        
	        if(printROC){
	        	for(WeightMatrix matrix : motifList){
		        	System.out.println("ROC:");
		        	finder = new MotifThresholdFinder(currgen, matrix, pos, neg, win);
			        finder.setROC(printROC);
			        double thres = finder.execute();
	        	}
	        }else{
		        System.out.println("Name\tMin\tMax\tBestThres\tBestF\tBestSn\tBestSp");
		        for(WeightMatrix matrix : motifList){
			        //Run the threshold finder
					//System.err.println("Initializing the threshold finder");
			    finder = new MotifThresholdFinder(currgen, matrix, pos, neg, win);
			        finder.setROC(printROC);
			       
			        //System.err.println("Finding the best threshold");
			        double bestThres = finder.execute();
			}
	        }
	       
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	//Constructors
    public MotifThresholdFinder(Genome g, WeightMatrix wm, ArrayList<Region> pos, ArrayList<Region> neg, int win){
		gen = g;
		posSet = pos;
		negSet = neg;
		motif=wm;
		window = win;
		if(wm==null){System.err.println("No motif specified");System.exit(1);}
		//Get the sequences
		SequenceGenerator seqgen = new SequenceGenerator();
		System.err.println("Fetching positive peak sequences");
		for(Region r : posSet){
		    posSeq.add(seqgen.execute(r));
		}
		System.err.println("Fetching negative peak sequences");
		for(Region r : negSet){
		    negSeq.add(seqgen.execute(r));
		}
	}
	
	public void setROC(boolean r){ROC = r;}
    public void setWin(int w){window=w;}

	//Find the motif-scoring threshold for the given specificity rate
	public double execute(){
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		double bestThres=0.0;
		
		//Find the scores for the positive sequences
		System.err.println("Scoring positive peaks with motif "+motif.getName());
		ArrayList<Double> posScore=new ArrayList<Double>();
		for(String s : posSeq){
			WeightMatrixScoreProfile profiler = scorer.execute(s);
			posScore.add(new Double(profiler.getMaxScore(profiler.getMaxIndex())));			
		}
		Collections.sort(posScore);

		//Find the scores for the negative sequences
		System.err.println("Scoring negative peaks with motif "+motif.getName());
		ArrayList<Double> negScore=new ArrayList<Double>();
		for(String s : negSeq){
			WeightMatrixScoreProfile profiler = scorer.execute(s);
			negScore.add(new Double(profiler.getMaxScore(profiler.getMaxIndex())));			
		}
		Collections.sort(negScore);
		
		
		//Find the score which maximizes performance
		bestThres=motif.getMaxScore();
		double bestF=0, bestSn=0, bestSp=0; int bestPerfIndex=0; 
		double posCount = (double)posScore.size(); double negCount = (double)negScore.size();
		
		if(ROC){
			System.out.println("i\tThreshold\tPerformance\tSn\tSp");
		}
		for(int d=posScore.size()-1; d>=0; d-=stepSize){
			double i = (double)(posScore.size()-d);
			double currThres = posScore.get(d);
			double currSn = (i/posCount);
			
			double negOverThres=0; 
			for(int e=negScore.size()-1; e>=0 && negScore.get(e)>=currThres; e-=stepSize){
				negOverThres+=stepSize;
			}
			double currSp = 1-(negOverThres/negCount);
			
			//Harmonic mean of Sn & Sp
			//double currF =2 * (currSn * currSp) / (currSn+currSp);

			//Maximum enrichment
			//double currF =(currSn) / (1-currSp);

			//PPV: TP/(TP+FP)
			//double currF=0;
			//if(currSp<1){
			//    double FP = (gen.getGenomeLength()/(double)window)*(1-currSp);
			//    currF = i/(i+FP);
			//}
			
			//F-beta
			double beta=1;
			double TP = i;
			double FP = (gen.getGenomeLength()/(double)window)*(1-currSp);
			double TN = (gen.getGenomeLength()/(double)window)*(currSp);
			double FN = posCount-i;
			double currF = ((1+beta*beta)*TP)/((1+beta*beta)*TP + (beta*beta)*FN + FP);

			//Accuracy
			//double currF = (TP+TN)/(TP+FP+TN+FN);

			if(currSp<1 && currF>bestF){
				bestSn=currSn;
				bestSp=currSp;
				bestPerfIndex=d;
				bestF=currF;
			}
			
			if(ROC){
			    System.out.println(i+"\t"+currThres+"\t"+currF+"\t"+currSn+"\t"+currSp);
			}
		}bestThres = posScore.get(bestPerfIndex);
		double max = motif.getMaxScore();
		double min = motif.getMinScore();
		System.out.println(motif.getName()+"\t"+min+"\t"+max+"\t"+bestThres+"\t"+bestF+"\t"+bestSn+"\t"+bestSp);
		
		return bestThres;
	}

	
	
	//Load a set of regions from a peak file
	public static ArrayList<Region> loadRegionsFromPeakFile(Genome gen,String filename, int win){
		ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line = reader.readLine(); //Ignore first line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>=3 && win!=-1){
	                PointParser pparser = new PointParser(gen);
	            	Point p = pparser.execute(words[2]);
	            	int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
                	int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
                	regs.add(r);
                }else if(words.length>=1){
	            	RegionParser parser = new RegionParser(gen);
	            	Region q = parser.execute(words[0]);
	            	if(win!=-1){
	            		int rstart = q.getMidpoint().getLocation()-(win/2)<1 ? 1:q.getMidpoint().getLocation()-(win/2);
	                	int rend = q.getMidpoint().getLocation()+(win/2)>gen.getChromLength(q.getChrom()) ? gen.getChromLength(q.getChrom()):q.getMidpoint().getLocation()+(win/2)-1;
	                	Region r = new Region(q.getGenome(), q.getChrom(), rstart, rend);
	                	if(r!=null){regs.add(r);}
	            	}else{
	            		if(q!=null){regs.add(q);}
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
	//Randomly pick a set of Regions
	public static ArrayList<Region> randomRegionPick(Genome gen, int numSamples, int sampleSize){
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
						validSamples++;
						regs.add(potential);
					}
				}total+=chromoSize[c];
			}				
		}
		return(regs);
	}
}
