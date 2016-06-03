package org.seqcode.projects.akshay.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.RepeatMaskedRegion;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.gse.utils.sequence.SequenceUtils;
import org.seqcode.projects.shaun.MotifAnalysisSandbox;


public class MotifToKmers {
	
	private GenomeConfig gcon;
	private List<Region> posSet;
	private List<Point> posPeaks=null;
	private List<Region> negSet;
	private List<Point> negPeaks = null;
	private int Kmin;
	private int Kmax;
	private WeightMatrix motif;
	//private double[] threshholds;
	private int threshlevel;
	private String posOutFile;
	private String negOutFile;
	
	
	
	
	
	public MotifToKmers (GenomeConfig gconf) {
		gcon = gcon;
	}
	
	
	//Mutators
	public void setPosSet(List<Region> pos){posSet = pos;}
	public void setNegSet(List<Region> neg){negSet = neg;}
	public void setKmin(int kmin){Kmin=kmin;}
	public void setKmax(int kmax){Kmax = kmax;}
	public void setMontif(WeightMatrix matrix){motif = matrix;}
	public void  setThreslevel(int t){threshlevel = t;}
	public void setPosOut(String po){posOutFile = po;}
	public void setNegOut(String no){negOutFile = no;}
	
	
	public static void main(String[] args) throws IOException, ParseException{
		ArgParser ap = new ArgParser(args);
		GenomeConfig gcon = new GenomeConfig(args);
		MotifToKmers runner = new MotifToKmers(gcon);
		boolean cache = ap.hasKey("seq");
		String seqpath="";
		if(cache){
			seqpath = ap.getKeyValue("seq");
		}
		int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
		String peaksFile = ap.getKeyValue("peaks");
		List<Point> peaksP = RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), ap.getKeyValue("peaksP"), win);
		List<Region> regsP = new ArrayList<Region>();
		for(Point p : peaksP){
			regsP.add(p.expand(win/2));
		}
		
		
		List<Point> peaksN = RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), ap.getKeyValue("peaksN"), win);
		List<Region> regsN = new ArrayList<Region>();
		for(Point p : peaksN){
			regsN.add(p.expand(win/2));
		}
		
		String outname = ap.getKeyValue("out");
		String posout = outname.concat("-pos.kmers");
		String negout = outname.concat("-neg.kmers");
		
		int thresLvl = ap.hasKey("threslevel") ? new Integer(ap.getKeyValue("threslevel")).intValue():10;
		
		int kmin = ap.hasKey("kmin") ? new Integer(ap.getKeyValue("kmin")).intValue():3;
		int kmax = ap.hasKey("kmax") ? new Integer(ap.getKeyValue("kmax")).intValue():8;
		
		String motiffile = ap.getKeyValue("motiffile");
		String backfile = ap.getKeyValue("back");
		
		List<WeightMatrix> matrixList = MotifAnalysisSandbox.loadMotifFromFile(motiffile, backfile, gcon.getGenome());
		WeightMatrix matrix = matrixList.get(0);
		
		runner.setPosSet(regsP);
		runner.setNegSet(regsN);
		runner.setKmin(kmin);
		runner.setKmax(kmax);
		runner.setMontif(matrix);
		runner.setThreslevel(thresLvl);
		runner.setPosOut(posout);
		runner.setNegOut(negout);
		
		runner.execute(cache, seqpath);
		
		
	}
	
	
	
	
	
	
	
	
	
	@SuppressWarnings("unchecked")
	public void execute(boolean useCache, String genPath) throws IOException{
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genPath);
		}
		
		FileWriter foutP = new FileWriter(posOutFile);
		FileWriter foutN = new FileWriter(negOutFile);
		StringBuilder headerSB = new StringBuilder();
		StringBuilder posSB = new StringBuilder();
		StringBuilder negSB = new StringBuilder();
		
		headerSB.append("Region");
		
		HashMap<String, List<Integer>> kmerCountsP = new HashMap<String,List<Integer>>();
		HashMap<String, List<Integer>> kmerCountsN = new HashMap<String,List<Integer>>();
		List<String> kmer_order=new ArrayList<String>();
		for(int k=Kmin; k<=Kmax; k++){ // iterating over different k-mer lengths
			int numK = (int)Math.pow(4, k);
			int[][] currKKmerCountsP = new int[posSet.size()][numK];
			int[][] currKKmerCountsN = new int[negSet.size()][numK];
			boolean[] isMotifKmer = new boolean[numK];
			
			List<Double> scores = new ArrayList<Double>();
			for(int i=0; i<numK; i++){ //iterating over k-mers of a given length and finding motifs that belong to a motif
				String currKmer = RegionFileUtilities.int2seq(i, k);
				String revCurrKmer = SequenceUtils.reverseComplement(currKmer);
				double scoref = scoreKmer(currKmer);
				double scoreR = scoreKmer(revCurrKmer);
				double score = scoref > scoreR ? scoref : scoreR;
				scores.add(score);
			}
			@SuppressWarnings("rawtypes")
			Comparator cmp = Collections.reverseOrder();  
			Collections.sort(scores,cmp);
			
			double threshold = scores.get((int)(threshlevel*scores.size()/100));
			
			for(int i=0; i<numK; i++){
				String currKmer = RegionFileUtilities.int2seq(i, k);
				String revCurrKmer = SequenceUtils.reverseComplement(currKmer);
				int revCurrKmerInd = RegionFileUtilities.seq2int(revCurrKmer);
				double scoref = scoreKmer(currKmer);
				double scoreR = scoreKmer(revCurrKmer);
				double score = scoref > scoreR ? scoref : scoreR;
				
				if(score >threshold && i < revCurrKmerInd  ){
					isMotifKmer[i] = true;
				}
				
			}
			
			for(int pr=0; pr<posSet.size(); pr++){
				for(int i=0; i <numK; i++){
					currKKmerCountsP[pr][i]=0;
				}
			}
			
			int prCounter=0;
			for(Region r : posSet){ // Count kmer frequencies at the positive set
				
				if(!kmerCountsP.containsKey(r.getLocationString())){
					kmerCountsP.put(r.getLocationString(),new ArrayList<Integer>());
				}
				String seq = seqgen.execute(r).toUpperCase();
				if(seq.contains("N"))
					continue;
				
				for(int i=0; i<(seq.length()-k+1); i++){
					String currK = seq.substring(i, i+k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					if(isMotifKmer[RegionFileUtilities.seq2int(currK)]){
						currKKmerCountsP[prCounter][RegionFileUtilities.seq2int(currK)]++;
					}
					if(isMotifKmer[RegionFileUtilities.seq2int(revCurrK)]){
						currKKmerCountsP[prCounter][RegionFileUtilities.seq2int(revCurrK)]++;
					}
				}
				prCounter++;
			}
			
			for(int nr=0; nr<negSet.size(); nr++){
				for(int i=0; i <numK; i++){
					currKKmerCountsN[nr][i]=0;
				}
			}
			
			int nrCounter = 0;
			for(Region r : negSet){
				if(!kmerCountsN.containsKey(r.getLocationString())){
					kmerCountsN.put(r.getLocationString(),new ArrayList<Integer>());
				}
				
				String seq = seqgen.execute(r).toUpperCase();
				if(seq.contains("N"))
					continue;
				
				for(int i=0; i<(seq.length()-k+1); i++){
					String currK = seq.substring(i, i+k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					if(isMotifKmer[RegionFileUtilities.seq2int(currK)]){
						currKKmerCountsN[nrCounter][RegionFileUtilities.seq2int(currK)]++;
					}
					if(isMotifKmer[RegionFileUtilities.seq2int(revCurrK)]){
						currKKmerCountsN[nrCounter][RegionFileUtilities.seq2int(revCurrK)]++;
					}
				}
				nrCounter++;
			}
			
			for(int i=0; i<numK; i++){ // Adding header
				if(isMotifKmer[i]){
					headerSB.append("\t");
					headerSB.append(RegionFileUtilities.int2seq(i, k));
				}
			}
			
			prCounter =0;
			for (Region r : posSet){
				for(int i=0; i<numK; i++){
					if(isMotifKmer[i]){
						kmerCountsP.get(r.getLocationString()).add(currKKmerCountsP[prCounter][i]);
						kmer_order.add(RegionFileUtilities.int2seq(i, k));
					}
				}
				prCounter++;
			}
			
			nrCounter = 0;
			for(Region r : negSet){
				for(int i=0; i<numK; i++){
					if(isMotifKmer[i]){
						kmerCountsN.get(r.getLocationString()).add(currKKmerCountsN[nrCounter][i]);
					}
				}
				nrCounter++;
			}
					
		}
		
		for(String rlocation : kmerCountsP.keySet()){
			posSB.append(rlocation);
			for(int c : kmerCountsP.get(rlocation)){
				posSB.append("\t");
				posSB.append(c);
			}
			posSB.append("\n");
		}
		
		for(String rlocation : kmerCountsN.keySet()){
			negSB.append(rlocation);
			for(int c : kmerCountsN.get(rlocation)){
				negSB.append("\t");
				negSB.append(c);
			}
			negSB.append("\n");
		}	
		
		foutP.write(headerSB.toString()+"\n");
		foutP.write(posSB.toString()+"\n");
		
		
		foutN.write(headerSB.toString()+"\n");
		foutN.write(negSB.toString()+"\n");

		foutP.close();
		foutN.close();
	}
	
	
	
	
	private double scoreKmer(String kmer){
		double maxScore=motif.getMinScore();
		char[] kmerChar = kmer.toCharArray();
		
		if(motif.length() < kmer.length()){
			maxScore = motif.getMinScore();
			return maxScore;
		}
		
		for(int i=0; i<motif.length()-kmer.length(); i++){
			float score = (float)0.0;
			for(int j=0; j< kmer.length(); j++){
				score += motif.matrix[i+j][kmerChar[j]];
			}
			
			if(maxScore <score){
				maxScore =score;
			}
		}
		
		return maxScore;
	}
	
	
	
	
	

}
