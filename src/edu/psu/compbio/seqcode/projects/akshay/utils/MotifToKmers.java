package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.RepeatMaskedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;

public class MotifToKmers {
	
	private GenomeConfig gcon;
	private List<Region> posSet;
	private List<Point> posPeaks=null;
	private List<Region> negSet;
	private List<Point> negPeaks = null;
	private int Kmin;
	private int Kmax;
	private WeightMatrix motif;
	private double[] threshholds;
	private String posOutFile;
	private String negOutFile;
	
	
	
	
	
	public MotifToKmers() {
		// TODO Auto-generated constructor stub
	}
	
	
	
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
		for(int k=Kmin; k<=Kmax; k++){ // iterating over different k-mer lengths
			int numK = (int)Math.pow(4, k);
			int[] kmerCounts = new int[numK];
			boolean[] isMotifKmer = new boolean[numK];
			
			for(int i=0; i<numK; i++){ //iterating over k-mers of a given length
				String currKmer = RegionFileUtilities.int2seq(i, k);
				String revCurrKmer = SequenceUtils.reverseComplement(currKmer);
				int revKmerIndex = RegionFileUtilities.seq2int(revCurrKmer);
				double scoref = scoreKmer(currKmer);
				double scoreR = scoreKmer(revCurrKmer);
				
				double score = scoref > scoreR ? scoref : scoreR;
				if(score >= threshholds[k-Kmin] && i < revKmerIndex){
					isMotifKmer[i] = true;
				}
			}
			
			for(Region r : posSet){
				posSB.append(r.getLocationString());
				for(int i=0; i<numK; i++){
					kmerCounts[i] = 0;
				}
				String seq = seqgen.execute(r).toUpperCase();
				if(seq.contains("N"))
					continue;
				
				for(int i=0; i<(seq.length()-k+1); i++){
					String currK = seq.substring(i, i+k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					if(isMotifKmer[RegionFileUtilities.seq2int(currK)]){
						kmerCounts[RegionFileUtilities.seq2int(currK)]++;
					}
					if(isMotifKmer[RegionFileUtilities.seq2int(revCurrK)]){
						kmerCounts[RegionFileUtilities.seq2int(revCurrK)]++;
					}
				}
			}
			
			for(int i=0; i<numK; i++){
				headerSB.append("\t");
				headerSB.append(RegionFileUtilities.int2seq(i, k));
				
			}
				//foutP.write("#"+getProgramName()+"\n");
			
			
			
			
			
			
			
			
			
		}
		
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
				score += motif.matrix[i][kmerChar[j]];
			}
			if(maxScore <score){
				maxScore =score;
			}
		}
		return maxScore;
	}
	
	
	
	
	

}
