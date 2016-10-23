package org.seqcode.motifs;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.motifdb.CountsBackgroundModel;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.sequence.RandomSequenceGenerator;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;


public class MarkovMotifThresholdFinder {
	private WeightMatrix motif = null;
	private MarkovBackgroundModel back;
	private ArrayList<String> seqSet = new ArrayList<String>(); 
	private static int numTest=100000;
	private int window=300;
	private boolean ROC=false;
	private boolean scored=false;
	private boolean seqGenerated=false;
	private ArrayList<Double> scores;
	
	public static void main(String[] args) throws IOException, ParseException {
		MarkovMotifThresholdFinder finder;
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||(!ap.hasKey("motifname")&&!ap.hasKey("motiffile"))) { 
            System.err.println("Usage:\n " +
                               "MarkovMotifThresholdFinder " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--motifname <weightmatrix name> "+
                               "--motifversion <weightmatrix version> " +
                               "--motiffile <file containing motifs> "+
                               "--back <background Markov model> "+
                               "--win <window of sequence around positive/negative points> "+
                               "--num <number of sequences to sample> " +
                               "--printroc ");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String motifversion=null;
        if(ap.hasKey("motifversion")){motifversion = ap.getKeyValue("motifversion");}
        String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        int numSim = 1000000;
        if(ap.hasKey("num")){
        	numSim = new Integer(ap.getKeyValue("num")).intValue();
        }
        boolean printROC= ap.hasKey("printroc");
        boolean loadFromFile = ap.hasKey("motiffile");
        
        
        try {
			//Load genome
			Species currorg = Species.getSpecies(species);
			//Genome currgen = currorg.getGenome(genome);

	        //Load the background model
	        MarkovBackgroundModel backMod;
	        if(backFile == null){
	          backMod = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(Genome.findGenome(genome)));
	        }else{
	        	backMod = BackgroundModelIO.parseMarkovBackgroundModel(backFile, Genome.findGenome(genome));
	        }
	        
	        //Pre-load the random sequences
	        ArrayList<String> randSeq = new ArrayList<String>();
	        RandomSequenceGenerator gen = new RandomSequenceGenerator(backMod);
			for(int i=0; i<numTest; i++){
				randSeq.add(gen.execute(win));
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
		        	finder = new MarkovMotifThresholdFinder(matrix, backMod, numSim);
			        if(win >0){finder.setWin(win);}
			        finder.setRandomSeq(randSeq);
			        finder.setROC(printROC);
			        double thres_1 = finder.execute(0.1);
	        	}
	        }else{
		        System.out.println("Name\tMin\tMax\tThres0.1\tThres0.05\tThres0.01\tThres0.005\tThres0.001");
		        for(WeightMatrix matrix : motifList){
			        //Run the threshold finder
					//System.err.println("Initializing the threshold finder");
			        finder = new MarkovMotifThresholdFinder(matrix, backMod, numSim);
			        if(win >0){finder.setWin(win);}
			        finder.setRandomSeq(randSeq);
			        finder.setROC(printROC);
			       
			        //System.err.println("Finding the best threshold");
			        double max = matrix.getMaxScore();
			        double min = matrix.getMinScore();
			        double thres_001 = finder.execute(0.001);
			        double thres_005 = finder.execute(0.005);
			        double thres_01 = finder.execute(0.01);
			        double thres_05 = finder.execute(0.05);
			        double thres_1 = finder.execute(0.1);
			        
			        System.out.println(matrix.getName()+"\t"+min+"\t"+max+"\t"+thres_1+"\t"+thres_05+"\t"+thres_01+"\t"+thres_005+"\t"+thres_001);
			        //System.out.println("Threshold for Sp=0.005:\t"+thres_005);
			        //System.out.println("Threshold for Sp=0.01:\t"+thres_01);
			        //System.out.println("Threshold for Sp=0.05:\t"+thres_05);
		        }
	        }
	       
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	//Constructors
	public MarkovMotifThresholdFinder(WeightMatrix wm, MarkovBackgroundModel markov){
		this(wm, markov, numTest);
	}
	public MarkovMotifThresholdFinder(WeightMatrix wm, MarkovBackgroundModel markov, int numSim){
		numTest=numSim;
		motif=wm;
		if(wm==null){System.err.println("No motif specified");System.exit(1);}
		back=markov;		
	}
	
	public void setNumTest(int nt){numTest=nt;}
	public void setWin(int w){window=w;}
	public void setROC(boolean r){ROC = r;}
	public void setRandomSeq(ArrayList<String> rand){if(rand.size()>0){seqSet = rand; seqGenerated=true;}}
	
	//Find the motif-scoring threshold for the given specificity rate
	public double execute(double Sp){
		if(Sp<0 || Sp>1){System.err.println("Invalid Sp value in MarkovMotifThreshold");System.exit(1);}

		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		double bestThres=0.0;
		
		//Find the scores for the random sequences
		if(!scored){
			//Generate the sequences first 
			//Simulate sequences using the markov background
			if(!seqGenerated){
				RandomSequenceGenerator gen = new RandomSequenceGenerator(back);
				for(int i=0; i<numTest; i++){
					seqSet.add(gen.execute(window));
				}
				seqGenerated=true;
			}
			
			scores=new ArrayList<Double>();
			for(String s : seqSet){
				WeightMatrixScoreProfile profiler = scorer.execute(s);
				scores.add(new Double(profiler.getMaxScore(profiler.getMaxIndex())));
			}
			Collections.sort(scores);
			scored=true;
		}
			
		//Find the score which corresponds to the required Specificity rate
		int index = (int)((double)scores.size()*(1-Sp));
		bestThres=scores.get(index);
		
		//Print an ROC if required
		if(ROC){
			System.out.println("i\tThreshold\tPerformance\tSp");
			int count=1;
			for(Double d : scores){
				double currThres = d.doubleValue();
				double currSp =(double)count/(double)scores.size(); 
				System.out.println(count+"\t"+currThres+"\t"+currSp);
				count++;
			}
		}
		return bestThres;
	}

	public Score2Sp getMotifROC(){
		ArrayList<Pair<Double,Double>> scoreVsSp = new ArrayList<Pair<Double,Double>>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		
		//Find the scores for the random sequences
		if(!scored){
			//Generate the sequences first 
			//Simulate sequences using the markov background
			if(!seqGenerated){
				RandomSequenceGenerator gen = new RandomSequenceGenerator(back);
				for(int i=0; i<numTest; i++){
					seqSet.add(gen.execute(window));
				}
				seqGenerated=true;
			}
			
			scores=new ArrayList<Double>();
			for(String s : seqSet){
				WeightMatrixScoreProfile profiler = scorer.execute(s);
				scores.add(new Double(profiler.getMaxScore(profiler.getMaxIndex())));
			}
			Collections.sort(scores);
			scored=true;
		}
			
		int count=1;
		for(Double d : scores){
			double currThres = d.doubleValue();
			double currSp =(double)count/(double)scores.size(); 
			scoreVsSp.add(new Pair<Double,Double>(currThres,currSp));
			count++;
		}
		return(new Score2Sp(scoreVsSp));
	}
}
