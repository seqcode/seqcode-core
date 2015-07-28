package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.seed.features.SuperEnrichedFeature;

public class KmerModelScanner {
	
	protected GenomeConfig gcon;
	protected double[] kmerweights;
	protected double[][] kmerpairwights;
	protected boolean isPairModel = false;
	protected List<Point> peaks;
	protected List<Region> regions;
	protected int k; // kmer length of the model
	protected int minM; // Minimum length to consider for motif finding
	protected int maxM; //Maximum length to consider for motif finding
	//protected int n; // Number of motifs to look at each 
	
	public KmerModelScanner(GenomeConfig gc) {
		gcon=gc;
	}
	
	
	
	// Settors
	public void setKmerWeights(double[] w){kmerweights=w;}
	public void setPeaks(List<Point>ps){peaks=ps;}
	public void setRegions(List<Region>rs){regions=rs;}
	public void setK(int K){k=K;}
	public void setminM(int m){minM=m;}
	public void setmaxM(int m){maxM=m;}
	public void setModelType(boolean isPair){isPair = isPairModel;}
	public void setKmerPairWeights(double[][] pw){kmerpairwights = pw;}
	
	public void printModelProfileMatrix(boolean useCache, String genpath,int slide){
		@SuppressWarnings("rawtypes")
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genpath);
		}
		for(Region r : regions){
			String seq = seqgen.execute(r).toUpperCase();
			if(seq.contains("N"))
				continue;
			StringBuilder sb = new StringBuilder();
			sb.append(r.getLocationString()+"\t");
			for(int i=0; i<(seq.length()-slide+1); i++){
				String subseq = seq.substring(i, i+slide);
				double score=0.0;
				for(int j=0; j<subseq.length()-k+1; j++){
					String currK = subseq.substring(j,j+k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					int currKInt = RegionFileUtilities.seq2int(currK);
					int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
					score = score+kmerweights[kmer];
				}			
				sb.append(Double.toString(score)+"\t");
			}
			sb.deleteCharAt(sb.length()-1);
			System.out.println(sb.toString());
		}
	}
	
	public void printKmerMountainComposition(boolean useCache, String genpath, double oddsThresh){
		@SuppressWarnings("rawtypes")
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genpath);
		}
		int[][] composition = new int[regions.size()][(int) Math.pow(4, (double)k)];
		String[] colnames = new String[(int) Math.pow(4, (double)k)];
		for(int c=0; c<colnames.length; c++){
			colnames[c] = RegionFileUtilities.int2seq(c, k);
		}
		String[] rownames = new String[regions.size()];
		int RegioinInd = 0;
		for(Region r : regions){
			String seq = seqgen.execute(r).toUpperCase();
			if(seq.contains("N")){
				rownames[RegioinInd] = r.getLocationString();
				for(int c =0;c<colnames.length; c++)
					composition[RegioinInd][c]=0;
				RegioinInd++;
				continue;
			}
			rownames[RegioinInd] = r.getLocationString();
			List<Pair<Region,Double>> mountains = new ArrayList<Pair<Region,Double>>();
			for(int l=minM; l<=maxM; l++){
				for(int i=0; i<(seq.length()-l+1); i++){
					String motif = seq.substring(i, i+l);
					double score=0.0;
					for(int j=0; j<motif.length()-k+1; j++){
						String currk = motif.substring(j, j+k);
						String revcurrk = SequenceUtils.reverseComplement(currk);
						int  currKInt = RegionFileUtilities.seq2int(currk);
						int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
						score = score+kmerweights[kmer];
					}
					Region hill = new Region(gcon.getGenome(),r.getChrom(),r.getStart()+i,r.getStart()+i+l-1);
					Iterator<Pair<Region,Double>> it = mountains.iterator();
					boolean add=true;
					while(it.hasNext() && add){
						Pair<Region,Double> pr = it.next();
						Region currHill = pr.car();
						Double currScore = pr.cdr();
						if(currHill.overlaps(hill) && currScore<score){
							it.remove();
							add=true;
						}else if(currHill.overlaps(hill) && currScore> score){
							add=false;
						}
					}
					if(add && score > oddsThresh){
						mountains.add(new Pair(hill,score));
					}
				}
				
			}
			
			for(Pair<Region,Double> hillpair: mountains){
				String motseq = seqgen.execute(hillpair.car()).toUpperCase();
				//int[] composition = new int[(int) Math.pow(4, (double)k)];
				//StringBuilder sb = new StringBuilder();
				//sb.append(hillpair.car().getLocationString()+"\t");
				for(int i=0; i<(motseq.length()-k+1);i++){
					String currK = motseq.substring(i, i+k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					int currKInt = RegionFileUtilities.seq2int(currK);
					int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt < revCurrKInt? currKInt:revCurrKInt;
					composition[RegioinInd][kmer]++;
				}
				//for(int i=0; i<(int) Math.pow(4, (double)k);i++){
				//	sb.append(composition[i]+"\t");
				//}
				//sb.deleteCharAt(sb.length()-1);
				//System.out.println(sb.toString());
			}
			RegioinInd++;
		}
	}
	
	public void printKmerMountains(boolean useCache, String genpath, double oddsThresh){
		@SuppressWarnings("rawtypes")
		SequenceGenerator seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genpath);
		}
		for(Region r : regions){
			String seq = seqgen.execute(r).toUpperCase();
			if(seq.contains("N"))
				continue;
			List<Pair<Region,Double>> mountains = new ArrayList<Pair<Region,Double>>();
			for(int l=minM; l<=maxM; l++){
				for(int i=0; i<(seq.length()-l+1); i++){
					String motif = seq.substring(i, i+l);
					double score=0.0;
					for(int j=0; j<motif.length()-k+1; j++){
						String currk = motif.substring(j, j+k);
						String revcurrk = SequenceUtils.reverseComplement(currk);
						int  currKInt = RegionFileUtilities.seq2int(currk);
						int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
						score = score+kmerweights[kmer];
					}
					if(isPairModel){
						for(int j=0; j<(motif.length()-k+1-1); j++){
							String currKP1 = motif.substring(j, j+k);
							String revCurrKP1 = SequenceUtils.reverseComplement(currKP1);
							int currKP1Int = RegionFileUtilities.seq2int(currKP1);
							int revCurrKP1Int = RegionFileUtilities.seq2int(revCurrKP1);
							int kmerP1 = currKP1Int < revCurrKP1Int? currKP1Int : revCurrKP1Int;
							for(int s=j+1;s<(motif.length()-k+1);s++){
								String currKP2 = motif.substring(s, s+k);
								String revCurrKP2 = SequenceUtils.reverseComplement(currKP2);
								int currKP2Int = RegionFileUtilities.seq2int(currKP2);
								int revCurrKP2Int = RegionFileUtilities.seq2int(revCurrKP2);
								int kmerP2 = currKP2Int < revCurrKP2Int ? currKP2Int: revCurrKP2Int;
								int x = kmerP1 < kmerP2 ? kmerP1: kmerP2;
								int y = kmerP1 < kmerP2 ? kmerP2: kmerP1;
								score = score+kmerpairwights[x][y];
							}
							
						}
					}
					
					Region hill = new Region(gcon.getGenome(),r.getChrom(),r.getStart()+i,r.getStart()+i+l-1);
					
					Iterator<Pair<Region,Double>> it = mountains.iterator();
					boolean add=true;
					while(it.hasNext() && add){
						Pair<Region,Double> pr = it.next();
						Region currHill = pr.car();
						Double currScore = pr.cdr();
						if(currHill.overlaps(hill) && currScore<score){
							it.remove();
							add=true;
						}else if(currHill.overlaps(hill) && currScore> score){
							add=false;
						}
					}
					if(add && score > oddsThresh){
						mountains.add(new Pair(hill,score));
					}
					
				}
			}
			
			for(Pair<Region,Double> hillpair: mountains){
				String motseq = seqgen.execute(hillpair.car()).toUpperCase();
				System.out.println(r.getLocationString()+"\t"+hillpair.car().getLocationString()+"\t"+hillpair.cdr()+"\t"+motseq+"\t"+SequenceUtils.reverseComplement(motseq));
			}	
		}
	}
	
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
			if(seq.contains("N"))
				continue;
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
		double[][] pws=null;
		boolean isPair=false;
		if(ap.hasKey("PairModel")){
			isPair = true;
		}
		BufferedReader reader = new BufferedReader(new FileReader(weightsFile));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\t");
            if(!isPair || !words[0].contains("-")){
            	int kind = RegionFileUtilities.seq2int(words[0]);
            	ws[kind] = Double.parseDouble(words[1]);
            }else if(isPair && words[0].contains("-")){
            	String[] subwords = words[0].split("-");
            	int kp1ind = RegionFileUtilities.seq2int(subwords[0]);
            	int kp2ind = RegionFileUtilities.seq2int(subwords[1]);
            	pws = new double[(int) Math.pow(4, k)][(int) Math.pow(4, k)];
            	pws[kp1ind][kp2ind] =  Double.parseDouble(words[1]);
            	pws[kp2ind][kp1ind] = Double.parseDouble(words[1]);
            }
        }
        boolean cache=false;
        String genPath = "";
        if(ap.hasKey("seq")){
        	cache = true;
        	genPath = ap.getKeyValue("seq");
        }
        
        KmerModelScanner scanner = new KmerModelScanner(gcon);
        scanner.setK(k);
        scanner.setKmerWeights(ws);
        if(isPair)
        	scanner.setKmerPairWeights(pws);
        scanner.setModelType(isPair);
        scanner.setmaxM(M);
        scanner.setminM(m);
        scanner.setPeaks(ps);
        scanner.setRegions(rs);
        if(ap.hasKey("scan")){
        	scanner.scanPeaks(cache, genPath);
        }
        if(ap.hasKey("printMountains")){
        	Double threshold = Double.parseDouble(ap.getKeyValue("oddsthresh"));
        	scanner.printKmerMountains(cache, genPath, threshold);
        }
        if(ap.hasKey("printProfileMatrix")){
        	int slidewin = Integer.parseInt(ap.getKeyValue("slide"));
        	scanner.printModelProfileMatrix(cache, genPath, slidewin);
        }
        if(ap.hasKey("printMountainsComposition")){
        	Double threshold = Double.parseDouble(ap.getKeyValue("oddsthresh"));
        	scanner.printKmerMountainComposition(cache, genPath, threshold);
        }
	}
	
	/**
	 * Removes all the the columns and rows that have all "0's" and prints the matrix 
	 * @param SparseMatrix
	 * @param colnames
	 * @param rownames
	 */
	public void printSimplifedMatrix(int[][] SparseMatrix, String[] colnames, String[] rownames){
		boolean[] columnstoRemove = new boolean[SparseMatrix[0].length];
		boolean[] rowstoRemove = new boolean[SparseMatrix.length];
		for(int c=0; c<columnstoRemove.length;c++){
			columnstoRemove[c] = true;
		}
		for(int r=0; r<rowstoRemove.length; r++){
			rowstoRemove[r] = true;
		}
		for(int r=0; r<rowstoRemove.length; r++){
			boolean RemoveThisRow=true;
			for(int c=0; c<columnstoRemove.length; c++){
				if(SparseMatrix[r][c] > 0){
					RemoveThisRow=false;
				}
			}
			rowstoRemove[r] = RemoveThisRow;
		}
		for(int c=0; c<columnstoRemove.length; c++){
			boolean RemoveThisCol=true;
			for(int r=0; r<rowstoRemove.length; r++){
				if(SparseMatrix[r][c] > 0){
					RemoveThisCol=false;
				}
			}
			columnstoRemove[c] = RemoveThisCol;
		}
		
		// Printing the header
		StringBuilder header = new StringBuilder();
		for(int c=0; c<columnstoRemove.length; c++){
			if(!columnstoRemove[c])
				header.append(colnames[c]+"\t");
		}
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		StringBuilder matrix = new StringBuilder();
		for(int r=0; r<rowstoRemove.length; r++){
			if(rowstoRemove[r])
				continue;
			matrix.append(rownames[r]+"\t");
			for(int c=0; c<columnstoRemove.length; c++){
				if(!columnstoRemove[c]){
					matrix.append(SparseMatrix[r][c]+"\t");
				}
			}
			matrix.deleteCharAt(matrix.length()-1);
			matrix.append("\n");
		}
		System.out.print(matrix.toString());
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
