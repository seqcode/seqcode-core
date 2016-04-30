package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.BinningParameters;
import edu.psu.compbio.seqcode.projects.akshay.utils.CGScoreProfile;
import edu.psu.compbio.seqcode.projects.akshay.utils.CGScorer;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

public class HillsAnalysisSandbox {
	
	// SeqUnwinder model parameters
	protected int Kmin =4;
	protected int Kmax =5;
	public HashMap<String,double[]> kmerweights = new HashMap<String,double[]>();
	/** Model names */
	protected List<String> kmerModelNames = new ArrayList<String>();
	public int numK;
	
	// Hills
	protected List<Region> hills = new ArrayList<Region>();
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	protected double motifScoreThresh = 0.2;
	
	//CG 
	protected double CGPercCutoff = 0.25;
	
	
	// SeqGenerator
	SequenceGenerator<Region> seqgen = null;
	
	//Settors
	public void setCGpercCutoff(double cg){CGPercCutoff = cg;}
	public void setMotifThresh(double motThresh){motifScoreThresh =motThresh;}
	public void setKmin(int kmin){Kmin = kmin;}
	public void setKmax(int kmax){Kmax = kmax;}
	public void setNumK(){
		numK = 0;
		for(int k=Kmin; k<=Kmax; k++ ){
			numK += (int)Math.pow(4, k);
		}
	}
	public void setHills(List<Region> modHills){hills = modHills;}
	// Load freq matrices
	public void loadMotifsFromFile(String filename, MarkovBackgroundModel b) {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifImport.setBackground(b);
		motifs.addAll(motifImport.readTransfacMatrices(filename));
	}
	public void setKmerWeights(String weightFileName) throws NumberFormatException, IOException{
		BufferedReader reader = new BufferedReader(new FileReader(weightFileName));
		String line;
		boolean header = false;
		while ((line = reader.readLine()) != null) {
			line = line.trim();
			String[] words = line.split("\t");
			if(words[0].charAt(0) == '#' || words[0].contains("Variable") || words[0].contains("Class")){
				header = true;
				for(int i=1; i<words.length; i++){
					kmerModelNames.add(words[i]);
					kmerweights.put(words[i], new double[numK]);
				}
			}else{
				if(!header){
					System.err.println("Please provide a header in the K-mer weight file");
					System.exit(1);
				}
				int ind = getKmerBaseInd(words[0]) + RegionFileUtilities.seq2int(words[0]);
				for(int i = 1; i < words.length; i++ ){
					kmerweights.get(kmerModelNames.get(i-1))[ind] = Double.parseDouble(words[i]);
				}
			}

		}reader.close();
	}

	// Gettors
	public int getKmerBaseInd(String kmer) {
		int baseInd = 0;
		for (int k = Kmin; k < kmer.length(); k++) {
			baseInd += (int) Math.pow(4, k);
		}
		return baseInd;
	}
	
	
	public HillsAnalysisSandbox(GenomeConfig gcon) {
		seqgen = gcon.getSequenceGenerator();
	}
	
	
	public void printHillStats(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_MotifScore\t"+modName+"_AA-KmerScore\t"+
				modName+"AT-KmerScore\t"+
				modName+"AC-KmerScore\t"+
				modName+"AG-KmerScore\t"+
				modName+"CA-KmerScore\t"+
				modName+"CC-KmerScore\t"+
				modName+"CG-KmerScore\t"+
				modName+"GA-KmerScore\t"+
				modName+"GC-KmerScore\t"+
				modName+"TA-KmerScore\t");

		//First get the k-mer for the given list of motifs
		HashSet<String> motifKmers = new HashSet<String>();
		for(WeightMatrix mot : motifs){
			motifKmers.addAll(WeightMatrix.getConsensusKmers(mot, Kmin, Kmax));
		}
		
		for(Region hillreg : hills){
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			double motifScore=0.0, AAScore=0.0, ATScore=0.0, ACScore =0.0, AGScore=0.0, CAScore=0.0, CCScore=0.0, CGScore=0.0, GAScore=0.0, GCScore=0.0, TAScore=0.0;
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					
					//Motif-Score
					if(motifKmers.contains(currk) || motifKmers.contains(revcurrk) || kmerweights.get(modName)[baseInd+kmer] < 0){
						motifScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AA-Score
					if(currk.contains("AA") || revcurrk.contains("AA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						AAScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AT-Score
					if(currk.contains("AT") || revcurrk.contains("AT") || kmerweights.get(modName)[baseInd+kmer] < 0){
						ATScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AC-Score
					if(currk.contains("AC") || revcurrk.contains("AC") || kmerweights.get(modName)[baseInd+kmer] < 0){
						ACScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AG-Score
					if(currk.contains("AG") || revcurrk.contains("AG") || kmerweights.get(modName)[baseInd+kmer] < 0){
						AGScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//CA-Score
					if(currk.contains("CA") || revcurrk.contains("CA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						CAScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//CC-Score
					if(currk.contains("CC") || revcurrk.contains("CC") || kmerweights.get(modName)[baseInd+kmer] < 0){
						CCScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//CG-Score
					if(currk.contains("CG") || revcurrk.contains("CG") || kmerweights.get(modName)[baseInd+kmer] < 0){
						CGScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//GA-Score
					if(currk.contains("GA") || revcurrk.contains("GA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						GAScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//GC-Score
					if(currk.contains("GC") || revcurrk.contains("GC") || kmerweights.get(modName)[baseInd+kmer] < 0){
						GCScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//TA-Score
					if(currk.contains("TA") || revcurrk.contains("TA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						TAScore += kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(motifScore)
					+"\t"+Double.toString(Math.round(AAScore*100)/100.0)+
					"\t"+Double.toString(Math.round(ATScore*100)/100.0)+
					"\t"+Double.toString(Math.round(ACScore*100)/100.0)+
					"\t"+Double.toString(Math.round(AGScore*100)/100.0)+
					"\t"+Double.toString(Math.round(CAScore*100)/100.0)+
					"\t"+Double.toString(Math.round(CCScore*100)/100.0)+
					"\t"+Double.toString(Math.round(CGScore*100)/100.0)+
					"\t"+Double.toString(Math.round(GAScore*100)/100.0)+
					"\t"+Double.toString(Math.round(GCScore*100)/100.0)+
					"\t"+Double.toString(Math.round(TAScore*100)/100.0));
		}

	}
	
	// Prints the CG hills that do not score highly for the given list of motifs (usually primary motifs)
	public void printCGhills(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_MotifScore\t"+"fractionCG");
		//First get the k-mer for the given list of motifs
		HashSet<String> motifKmers = new HashSet<String>();
		for(WeightMatrix mot : motifs){
			motifKmers.addAll(WeightMatrix.getConsensusKmers(mot, Kmin, Kmax));
		}
		
		for(Region hillreg : hills){
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			double score=0.0;
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					if(motifKmers.contains(currk) || motifKmers.contains(revcurrk)){
						int  currKInt = RegionFileUtilities.seq2int(currk);
						int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
						int baseInd = this.getKmerBaseInd(currk);
						score = score+kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			//Now if th motif score is less check for CGs
			if(score < motifScoreThresh){
				BinningParameters params = new BinningParameters(hillreg.getWidth(), hillreg.getWidth());
				CGScorer scorer = new CGScorer(params,seqgen);
				int MaxCG = (int)hillreg.getWidth()/2;
				CGScoreProfile profile = scorer.execute(hillreg);
				int total = profile.getCGcount();
				double perc = total/(double)MaxCG;
				if(perc > CGPercCutoff)
					System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score)+"\t"+Double.toString(perc));
			}
			
		}
	
	}
	
	public void scoreHillsWithCGKmers(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_CGKmerScore");
		
		for(Region hillreg : hills){
			double score=0.0;
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					
					if(currk.contains("CG") || revcurrk.contains("CG") || kmerweights.get(modName)[baseInd+kmer]<0){
						score = score+kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score));
		}
		
		
	}
	
	
	public void scoreHills(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_MotifScore");
		
		for(Region hillreg : hills){
			double score=0.0;
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					score += kmerweights.get(modName)[baseInd+kmer];
				}
			}
			System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score));
		}
	}
		

	
	// Prints all the hills that are high scoring for the given list of motifs
	public void printMotifHills(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_Score\t");
		//First get the k-mer for the given list of motifs
		HashSet<String> motifKmers = new HashSet<String>();
		for(WeightMatrix mot : motifs){
			motifKmers.addAll(WeightMatrix.getConsensusKmers(mot, Kmin, Kmax));
		}

		for(Region hillreg : hills){
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			double score=0.0;
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					if(motifKmers.contains(currk) || motifKmers.contains(revcurrk) || kmerweights.get(modName)[baseInd+kmer] < 0){
						score = score+kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			
			if(score > motifScoreThresh){
				System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score));
			}
		}
	}
	
	
	public static void main(String[] args) throws IOException, ParseException{
		ArgParser ap = new ArgParser(args);
		GenomeConfig gcon = new GenomeConfig(args);
		HillsAnalysisSandbox runner = new HillsAnalysisSandbox(gcon);

		String backFile =ap.hasKey("back") ? ap.getKeyValue("back"):null;

		MarkovBackgroundModel back;

		if(backFile == null){
			back = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gcon.getGenome()));
		}else{
			back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gcon.getGenome());
		}

		// Length of the smallest K-mer in the K-mer models
		int minK = Args.parseInteger(args, "Kmin", 4);
		runner.setKmin(minK);

		// Length of the largest K-mer in the K-mer models
		int maxK = Args.parseInteger(args, "Kmax", 5);
		runner.setKmax(maxK);

		runner.setNumK();

		// K-mer models file / weights file
		String weights = Args.parseString(args, "weights", null);
		if (weights == null) {
			System.err.println("Provide weights file");
			System.exit(1);
		}
		runner.setKmerWeights(weights);
		
		// Load hills
		String hillsFile = ap.getKeyValue("hillsFile");
		if(hillsFile == null){
			System.err.println("Provide hills file");
			System.exit(1);
		}
		List<Region> regs = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), hillsFile, -1);
		runner.setHills(regs);

		// Now load K-mer models
		String motifFile = ap.getKeyValue("motiffile");
		if(ap.hasKey("motiffile"))
			runner.loadMotifsFromFile(motifFile,back);
		
		double cgThresh = Args.parseDouble(args, "cgThresh", 0.25);
		runner.setCGpercCutoff(cgThresh);
		
		double motThresh = Args.parseDouble(args, "motThresh", 0.2);
		runner.setMotifThresh(motThresh);
		
		String modName = ap.getKeyValue("modname");
		if(modName == null){
			System.err.println("Which modelname to scan? Please provide and re-run");
			System.exit(1);
		}
			
		if(ap.hasKey("scoreHills"))
			runner.scoreHills(modName);
		if(ap.hasKey("printCGhills"))
			runner.printCGhills(modName);
		if(ap.hasKey("printMotHills"))
			runner.printMotifHills(modName);
		if(ap.hasKey("scoreHillsWithCGKmers"))
			runner.scoreHillsWithCGKmers(modName);
		if(ap.hasKey("stats"))
			runner.printHillStats(modName);

	}

}
