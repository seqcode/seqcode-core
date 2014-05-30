package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredStrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;


import edu.psu.compbio.seqcode.gse.utils.probability.Hypergeometric;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.shaun.MotifAnalysisSandbox;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

public class KmerMapper {
	
	public List<StrandedPoint> points = new ArrayList<StrandedPoint>();
	public List<StrandedRegion> regions = new ArrayList<StrandedRegion>();
	public int k;
	public SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	public int[][][] MapMatrix;
	public int[][][] PairWiseMatrix;
	public int winSize;
	public Genome gen;
	public WeightMatrix motif;
	public WeightMatrixScorer scorer;
	public double motifThres;
	public int[][][] PairWiseDistMatrix;
	
	
	
	public List<String> motif_kmers = new ArrayList<String>();
	
	public KmerMapper(Genome g, int win, double threshold, WeightMatrix mot, int ksize) {
		this.gen=g;
		this.winSize=win;
		this.motifThres = threshold;
		this.motif = mot;
		this.k = ksize;
	}
	
	public KmerMapper(Genome g, int win, int ksize){
		this.gen = g;
		this.winSize = win;
		this.k = ksize;
		this.motif = null;
		this.motifThres = Double.MIN_VALUE;
	}
	
	//Settors
	
	public void setMotifRegions(String peaksFileName,String SeqPathFile){
		List<Point> tempPoints = Utilities.loadPeaksFromPeakFile(gen, peaksFileName, winSize);
		List<Region> tempRegions = Utilities.loadRegionsFromPeakFile(gen, peaksFileName, winSize);
	
		seqgen.useCache(true);
		seqgen.setGenomePath(SeqPathFile);
		scorer = new WeightMatrixScorer(motif,seqgen);
			
		for(int i=0 ; i<tempRegions.size(); i++){
			Point p= tempPoints.get(i);
			Region r = tempRegions.get(i);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			int bestMotifIndex = profiler.getMaxIndex();
			double bestMotifScore = profiler.getMaxScore(bestMotifIndex);
			if(bestMotifScore >= motifThres){
				int closestDist = Integer.MAX_VALUE;
				int closestIndex = -1;
				char closestStrand = '\0';
				for(int z=0; z<r.getWidth()-motif.length()+1; z++){
					double currScore= profiler.getMaxScore(z);
					if(currScore>=motifThres){
						int motifCenterCoord = z+(motif.length()/2)+r.getStart();
						int dist = Math.abs(p.getLocation() - motifCenterCoord);
						if(dist<closestDist){
							closestDist = dist;
							closestIndex = z;
							closestStrand = profiler.getMaxStrand(z);
						}
					}
				}
				
				if(closestStrand == '+'){
					Point addP = new Point(gen,r.getChrom(),r.getStart()+closestIndex);
					points.add(new StrandedPoint(addP,closestStrand));
				}else{
					Point addP = new Point(gen,r.getChrom(),r.getStart()+closestIndex+motif.length()-1);
					points.add(new StrandedPoint(addP,closestStrand));
				}
			}
		}
		for(StrandedPoint sp : points){
			regions.add(new StrandedRegion(sp.expand(winSize/2),sp.getStrand()));
		}
	}
	
	public void setRegions(String peaksFileName,String SeqPathFile){
		List<Point> tempPoints = Utilities.loadPeaksFromPeakFile(gen, peaksFileName, winSize);
		List<Region> tempRegions = Utilities.loadRegionsFromPeakFile(gen, peaksFileName, winSize);
		for(Point p : tempPoints){
			points.add(new StrandedPoint(p,'+'));
		}
		for(Region r : tempRegions){
			regions.add(new StrandedRegion(r,'+'));
		}
	}
	
	public void setMapMatrix(String SeqPathFile){
		int numk = (int)Math.pow(4, k);
		this.MapMatrix = new int[numk][regions.size()][winSize];
		seqgen.useCache(true);
		seqgen.setGenomePath(SeqPathFile);
		
		for(int i=0; i<regions.size(); i++){
			String seq = seqgen.execute(regions.get(i)).toUpperCase();
			if(seq.contains("N"))
				continue;
			String seqOriented = (regions.get(i).getStrand() == '+') ? seq : SequenceUtils.reverseComplement(seq);
			for(int j=0; j<(seq.length()-k+1); j++){
				String currK = seqOriented.substring(j, j+k);
				String revCurrK = SequenceUtils.reverseComplement(currK);
				int currKInt = Utilities.seq2int(currK);
				int revCurrKInt = Utilities.seq2int(revCurrK);
				int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
				MapMatrix[kmer][i][j] = 1;
			}
		}
		
	}
	
	public void setPairWiseMatrix(String SeqPathFile){
		int numk = (int)Math.pow(4, k);
		this.PairWiseMatrix = new int[numk][regions.size()][winSize];
		seqgen.useCache(true);
		seqgen.setGenomePath(SeqPathFile);
		for(int i=0; i<regions.size(); i++){
			String seq = seqgen.execute(regions.get(i)).toUpperCase();
			if(seq.contains("N"))
				continue;
			String seqOriented = (regions.get(i).getStrand() == '+') ? seq : SequenceUtils.reverseComplement(seq);
			for(int j=0; j<(seq.length()-k+1); j++){
				String currK = seqOriented.substring(j, j+k);
				String revCurrK = SequenceUtils.reverseComplement(currK);
				int currKInt = Utilities.seq2int(currK);
				int revCurrKInt = Utilities.seq2int(revCurrK);
				int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
				PairWiseMatrix[kmer][i][j] = kmer;
			}
		}
		
	}
	
	public void setMotifKmers(int percCutoff){
		for(int i=0; i< MapMatrix.length; i++){
			String kmer = Utilities.int2seq(i, k);
			String revkmer = SequenceUtils.reverseComplement(kmer);
			int[][] kmerMap = MapMatrix[i];
			int revi = Utilities.seq2int(revkmer);
			if(i>revi)
				continue;
			boolean isMotifKmer = false;
			int midP = winSize/2;
			for(int j=midP; j<midP+motif.length()-k; j++){
				int rowSum = 0;
				for(int h=0; h< kmerMap.length; h++){
					rowSum = rowSum +kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff)
					isMotifKmer = true;
			}
			if(isMotifKmer)
				motif_kmers.add(kmer);
		}
	}
	
	public void setPairWiseDistMatrix(List<String> kmerSet){
		int numk = (int)Math.pow(4, k);
		this.PairWiseDistMatrix =  new int[this.winSize][numk][numk];
		for(int i=0; i<kmerSet.size(); i++){
			for(int j=i; j<kmerSet.size(); j++){
				int kmerXind = (Utilities.seq2int(kmerSet.get(i)) > Utilities.seq2int(kmerSet.get(j))) ? Utilities.seq2int(kmerSet.get(i)) :  Utilities.seq2int(kmerSet.get(j));
				int kmerYind =  (Utilities.seq2int(kmerSet.get(i)) > Utilities.seq2int(kmerSet.get(j))) ? Utilities.seq2int(kmerSet.get(j)) : Utilities.seq2int(kmerSet.get(i));
				for(int l=0; l<regions.size(); i++){ // over all locations 
					int[] Xmap = this.PairWiseMatrix[kmerXind][l];
					int[] Ymap = this.PairWiseMatrix[kmerYind][l];
					int[] XYmap = new int[Xmap.length];
					for(int k=0; k<XYmap.length; k++){
						XYmap[k] = Xmap[k]+Ymap[k];
					}
					int distance = 0;
					boolean counting = false;
					int leftInd = 0;
					for(int k=0; k<XYmap.length; k++){
						if(!counting && XYmap[k] != 0){
							counting = true;
							leftInd = XYmap[k];
							distance = 1;
						}
						if(counting && XYmap[k] == 0){
							distance++;
						}
						if(counting && XYmap[k] !=0 && XYmap[k] == leftInd){
							if(kmerXind == kmerYind){
								this.PairWiseDistMatrix[distance][kmerXind][kmerYind]++;
								distance = 1;
								leftInd = XYmap[k];
							}
							else{
								distance =1;
								leftInd = XYmap[k];
							}
						}
						if(counting && XYmap[k] !=0 && XYmap[k] != leftInd){
							this.PairWiseDistMatrix[distance][kmerXind][kmerYind]++;
							distance = 1;
							leftInd = XYmap[k];
						}
					}
				}
			
			}
		}
		
	}
	
	
	// Calculators
	public void printMotifKmerPvalue(List<String> kmerSet){
		int popSize = 136;
		int motSize = motif_kmers.size();
		int markedSize = 0;
		for(int i=0; i<kmerSet.size(); i++){
			for(int j=0; j< motif_kmers.size(); j++){
				if(kmerSet.get(i).equals(motif_kmers.get(j))){
					markedSize++;
					continue;
				}
			}
		}
		
		Hypergeometric tester = new Hypergeometric();
		double log_pvalue = tester.log_hypgeomPValue(popSize, motSize, kmerSet.size(), markedSize);
		double prob = Math.exp(log_pvalue);
		System.out.println(popSize+" "+motSize+" "+kmerSet.size()+" "+markedSize);
		System.out.println("Prob: "+prob );

	}

	
		
	public void printInformativeKmersFromSet(List<String> kmerSet, int percCutoff){
		for(int i=0; i<kmerSet.size(); i++){
			int kmerID = Utilities.seq2int(kmerSet.get(i));
			int[][] kmerMap = MapMatrix[kmerID];
			boolean pass = false;
			for(int j=0; j<kmerMap[0].length; j++){ // over all positions
				int rowSum =0;
				for(int h=0; h<kmerMap.length; h++){ // over all peaks
					rowSum = rowSum + kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff){
					pass = true;
				}
			}
			if(pass){
				System.out.println(kmerSet.get(i));
			}
		}
	}
	
	public void printInformativeKmers(int percCutoff){
		for(int i=0; i<MapMatrix.length; i++){
			int[][] kmerMap = MapMatrix[i];
			boolean pass = false;
			for(int j=0; j<kmerMap[0].length; j++){
				int rowSum =0;
				for(int h=0; h<kmerMap.length; h++){
					rowSum = rowSum + kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff){
					pass = true;
				}
			}
			if(pass){
				System.out.println(Utilities.int2seq(i, k));
			}
			
		}
	}
	
	public void printInformativeKmerWithinDistance(int percCutoff, int distance){
		for(int i=0; i<MapMatrix.length; i++){
			int[][] kmerMap = MapMatrix[i];
			boolean pass = false;
			int midP = (int)(kmerMap[0].length/2);
			for(int j=(midP-distance); j<(midP+distance); j++){ // over all locations
				int rowSum =0;
				for(int h=0; h<kmerMap.length; h++){ // over all peaks
					rowSum = rowSum + kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff){
					pass = true;
				}
			}
			if(pass){
				System.out.println(Utilities.int2seq(i, k));
			}
		}
		
	}
	
	public void printInformativeKmerWithinDistanceFromSet(List<String> kmerSet, int percCutoff, int distance){
		for(int i=0; i<kmerSet.size(); i++){
			int kmerID = Utilities.seq2int(kmerSet.get(i));
			int[][] kmerMap = MapMatrix[kmerID];
			boolean pass = false;
			int midP = (int)(kmerMap[0].length/2);
			for(int j=(midP-distance); j<(midP+distance); j++){ // over all locations
				int rowSum =0;
				for(int h=0; h<kmerMap.length; h++){ // over all peaks
					rowSum = rowSum + kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff){
					pass = true;
				}
			}
			if(pass){
				System.out.println(kmerSet.get(i));
			}
		}
	}
	
	/**
	 * excludes k-mers coming from motifs that have been centered
	 * @param percCutoff
	 * @param distance
	 */
	public void printInformativeKmerFlanks(int percCutoff, int distance){
		for(int i=0; i<MapMatrix.length; i++){
			int[][] kmerMap = MapMatrix[i];
			boolean pass = false;
			int midP = (int)(kmerMap[0].length/2);
			for(int j=(midP-distance); j<(midP+distance); j++){ // over all locations
				int rowSum =0;
				for(int h=0; h<kmerMap.length; h++){ // over all peaks
					rowSum = rowSum + kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff && (j<midP || j> (midP+motif.length()-k))){
					pass = true;
				}
			}
			if(pass){
				System.out.println(Utilities.int2seq(i, k));
			}
		}
	}
	
	public void printInformativeKmerFlanksFromSet(List<String> kmerSet, int percCutoff, int distance){
		for(int i=0; i<kmerSet.size(); i++){
			int kmerID = Utilities.seq2int(kmerSet.get(i));
			int[][] kmerMap = MapMatrix[kmerID];
			boolean pass = false;
			int midP = (int)(kmerMap[0].length/2);
			for(int j=(midP-distance); j<(midP+distance); j++){ // over all locations
				int rowSum =0;
				for(int h=0; h<kmerMap.length; h++){ // over all peaks
					rowSum = rowSum + kmerMap[h][j];
				}
				if((int)((rowSum*100)/kmerMap.length) > percCutoff && (j<midP || j> (midP+motif.length()-k))){
					pass = true;
				}
			}
			if(pass){
				System.out.println(kmerSet.get(i));
			}
		}
	}
	
	public static void main(String[] args) throws IOException, NotFoundException, ParseException{
		ArgParser ap = new ArgParser(args);
		System.err.println("Usage:\n"
		);
		
		String kmerFile = ap.getKeyValue("kmers");
		File kmerF = new File(kmerFile);
		BufferedReader reader = new BufferedReader(new FileReader(kmerF));
		String line;
		List<String> kmers = new ArrayList<String>();
		while((line = reader.readLine()) != null){
			line = line.trim();
			line.replace("\n", "");
			kmers.add(line);
		}
		reader.close();
		
		boolean printInfKmersS = false, printInfKmers = false, printInfKmersWDis = false, printInfKmersWDisS = false,
				printInfKmerFlanks = false, printInfKmerFlanksS = false;
		
		printInfKmersS = ap.hasKey("printInfKmersS");
		printInfKmers = ap.hasKey("printInfKmers");
		printInfKmersWDis = ap.hasKey("printKmersInfWDis");
		printInfKmersWDisS = ap.hasKey("printInfKmersWDisS");
		printInfKmerFlanks = ap.hasKey("printInfFlanks");
		printInfKmerFlanksS = ap.hasKey("printInfFlanksS");
		boolean printMotifKmerPvalue = ap.hasKey("motifKmerPval");
		
		
		String SeqFilePath = ap.getKeyValue("seq");
		int winSize = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue() : 200;
		Genome gen = null;
		if(ap.hasKey("species")){
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			if(pair != null){
				
				gen = pair.cdr();
			}
		}else{
			if(ap.hasKey("geninfo") || ap.hasKey("g")){
				//Make fake genome... chr lengths provided
				String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
				gen = new Genome("Genome", new File(fName), true);
			}else{
				gen = null;
			}
		}
		
		String motiffile, backfile;
		
		motiffile = ap.getKeyValue("motiffile");
		backfile = ap.getKeyValue("back");
		
		WeightMatrix matrix = null;
		matrix = MotifAnalysisSandbox.loadMotifFromFile(motiffile, backfile, gen).get(0);
		double minscore = (ap.hasKey("minscore")) ? new Double(ap.getKeyValue("minscore")).doubleValue() : 6.4;
		int ksize = ap.hasKey("ksize") ? new Integer(ap.getKeyValue("ksize")).intValue() : 4;
		String peaksFile = ap.getKeyValue("peaks");
		
		KmerMapper mapper = new KmerMapper(gen,winSize,minscore,matrix,ksize);
		mapper.setMotifRegions(peaksFile, SeqFilePath);
		mapper.setMapMatrix(SeqFilePath);
		
		
		int percCutoff = ap.hasKey("percentAlign") ? new Integer(ap.getKeyValue("percentAlign")).intValue() : 30;
		int distance = ap.hasKey("flankSize")? new Integer(ap.getKeyValue("flankSize")).intValue() : 30;
		
		mapper.setMotifKmers(percCutoff);
		
		if(printMotifKmerPvalue){
			System.out.println("HyperGeometric test result: ");
			mapper.printMotifKmerPvalue(kmers);
			
		}
		if(printInfKmers){
			System.out.println("Printing all Kmers that are aligned with respect to the primary motifs: ");
			mapper.printInformativeKmers(percCutoff);
		}
		
		if(printInfKmersS){
			System.out.println("Printing Kmers that are aligned with respect to the primary motif and have high weight in SVM classification");
			mapper.printInformativeKmersFromSet(kmers, percCutoff);
		}
		
		if(printInfKmersWDis){
			System.out.println("Printing all Kmers that are aligned with respect to the primary motifs and are within a small distance from the motif");
			mapper.printInformativeKmerWithinDistance(percCutoff, distance);
		}
		
		if(printInfKmersWDisS){
			System.out.println("Printing all Kmers that are aligned with respect to the primary motifs and are within a small distance from the motif and also have high SVM weight");
			mapper.printInformativeKmerWithinDistanceFromSet(kmers, percCutoff, distance);
		}
		if(printInfKmerFlanks){
			 System.out.println("Printing informative kmer flanks that don't come from the primary motif");
			 mapper.printInformativeKmerFlanks(percCutoff, distance);
		}
		
		if(printInfKmerFlanksS){
			System.out.println("Printing informative kmer flanks that don't come from the primart motif and have a high SVM score");
			mapper.printInformativeKmerFlanksFromSet(kmers, percCutoff, distance);
		}
		
	}

}