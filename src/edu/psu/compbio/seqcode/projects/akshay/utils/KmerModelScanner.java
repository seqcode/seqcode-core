package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class KmerModelScanner {
	
	
	protected double[] kmerweights;
	protected List<Point> peaks;
	protected List<Region> regions;
	protected int k; // kmer length of the model
	protected int minM; // Minimum length to consider for motif finding
	protected int maxM; //Maximum length to consider for motif finding
	//protected int n; // Number of motifs to look at each 
	
	public KmerModelScanner() {
		
	}
	
	
	
	// Settors
	public void setKmerWeights(double[] w){kmerweights=w;}
	public void setPeaks(List<Point>ps){peaks=ps;}
	public void setRegions(List<Region>rs){regions=rs;}
	public void setK(int K){k=K;}
	public void setminM(int m){minM=m;}
	public void setmaxM(int m){maxM=m;}
	
	public void scanPeaks(boolean useCache, String genpath){
		@SuppressWarnings("rawtypes")
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genpath);
		}
		
		for(Region r: regions){
			String seq = seqgen.execute(r).toUpperCase();
			double maxScore = Double.MIN_VALUE;
			String bestMotif = "";
			for(int l=minM; l<=maxM; l++){
				for(int i=0; i<(seq.length()-l+1);i++){
					String motif = seq.substring(i, i+l);
					double score =0.0;
					for(int j=0;j<motif.length()-k+1; j++){
						String currk =motif.substring(j, j+k);
						String revcurrk = SequenceUtils.reverseComplement(currk);
						int  currKInt = RegionFileUtilities.seq2int(currk);
						int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
						score = score+kmerweights[kmer];
					}
					if(score>maxScore){
						maxScore = score;
						bestMotif = motif;
					}
				}
			}
			System.out.println(r.getLocationString()+"\t"+bestMotif+"\t"+SequenceUtils.reverseComplement(bestMotif)+Double.toString(maxScore));
		}
		
		seqgen.clearCache();
		
	}
	
	
	public static void main(String[] args) throws IOException{
		GenomeConfig gcon = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		
		int k = Args.parseInteger(args, "k", 4);
		int m = Args.parseInteger(args, "m", 5);
		int M = Args.parseInteger(args, "M", 10);
		String peaksFile = ap.getKeyValue("peaks");
		int win = Args.parseInteger(args, "win", 150);
		List<Point> ps = RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), peaksFile, win);
		List<Region> rs = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), peaksFile, win);
		String weightsFile = Args.parseString(args, "weights", "");
		double[] ws = new double[(int) Math.pow(4, k)];
		BufferedReader reader = new BufferedReader(new FileReader(weightsFile));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\t");
            int kind = RegionFileUtilities.seq2int(words[0]);
            ws[kind] = Double.parseDouble(words[1]);
        }
        boolean cache=false;
        String genPath = "";
        if(ap.hasKey("seq")){
        	cache = true;
        	genPath = ap.getKeyValue("seq");
        }
        
        KmerModelScanner scanner = new KmerModelScanner();
        scanner.setK(k);
        scanner.setKmerWeights(ws);
        scanner.setmaxM(M);
        scanner.setminM(m);
        scanner.setPeaks(ps);
        scanner.setRegions(rs);
        
        scanner.scanPeaks(cache, genPath);
            
		
	}
	
	
	
	
	public int compareKmerModels(int[][] a, int[][] b) throws IncorrectComparision{
		int score=0;
		if(a.length != b.length){
			throw new IncorrectComparision("Unequal model lenghts");
		}else{
			for(int r=0; r<a.length; r++){
				for(int c=0; c<=a.length; c++){
					if(a[r][c] != b[r][c])
					score = score++;
				}
			}
		}		
		return score;
		
	}
	
	public class IncorrectComparision extends Exception{
		IncorrectComparision(String s){
			super(s);
		}
	}
	
	
}
