package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredStrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

public class KmerMapper {
	
	public List<Point> points;
	public List<Region> regions;
	public int k;
	public SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	public int[][][] MapMatrix;
	public int winSize;
	public Genome gen;
	public WeightMatrix motif;
	public WeightMatrixScorer scorer;
	public int motifThres;
	
	public KmerMapper(Genome g, int win, int threshold) {
		this.gen=g;
		this.winSize=win;
		this.motifThres = threshold;
	}
	
	//Settors
	
	public void setRegions(boolean shift, String peaksFileName,String SeqPathFile){
		List<Point> tempPoints = Utilities.loadPeaksFromPeakFile(gen, peaksFileName, winSize);
		List<Region> tempRegions = Utilities.loadRegionsFromPeakFile(gen, peaksFileName, winSize);
		if(shift){
			
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
					
					for(int z=0; z<r.getWidth()-motif.length()+1; z++){
						double currScore= profiler.getMaxScore(z);
						if(currScore>=motifThres){
							int motifCenterCoord = z+(motif.length()/2)+r.getStart();
							int dist = Math.abs(p.getLocation() - motifCenterCoord);
							if(dist<closestDist){
								closestDist = dist;
								closestIndex = z;
							}
						}
					}
					
					Region hitreg = new Region(gen,r.getChrom(), r.getStart()+closestIndex, r.getStart()+closestIndex+motif.length()-1);
					points.add(new Point(gen,hitreg.getChrom(),hitreg.getStart()));
				}
			}
			for(Point p : points){
				regions.add(p.expand(winSize/2));
			}
			
		}else{
			points = tempPoints;
			regions = tempRegions;
		}
		
	}
	
	
	public void setMapMatrix(String SeqPathFile){
		int numk = (int)Math.pow(4, k);
		this.MapMatrix = new int[numk][points.size()][winSize];
		seqgen.useCache(true);
		seqgen.setGenomePath(SeqPathFile);
		
		for(int i=0; i<regions.size(); i++){
			String seq = seqgen.execute(regions.get(i)).toUpperCase();
			if(seq.contains("N"))
				continue;
			for(int j=0; j<(seq.length()-k+1); j++){
				String currK = seq.substring(j, j+k);
				String revCurrK = SequenceUtils.reverseComplement(currK);
				int currKInt = Utilities.seq2int(currK);
				int revCurrKInt = Utilities.seq2int(revCurrK);
				int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
				boolean reverseComplement = currKInt > revCurrKInt;
				if(reverseComplement)
					MapMatrix[kmer][i][j+k-1] = 1;
				else
					MapMatrix[kmer][i][j]=1;
			}
		}
		
	}
	
	
	// Calculators
	
	public void printInformativeKmersFromSet(String[] kmerSet, int percCutoff){
		for(int i=0; i<kmerSet.length; i++){
			int kmerID = Utilities.seq2int(kmerSet[i]);
			int[][] kmerMap = MapMatrix[kmerID];
			boolean pass = false;
			for(int j=0; j<kmerMap[0].length; j++){
				int colSum =0;
				for(int k=0; k<kmerMap.length; k++){
					colSum = colSum + kmerMap[j][k];
				}
				if((int)((colSum*100)/kmerMap.length) > percCutoff){
					pass = true;
				}
			}
			if(pass){
				System.out.println(kmerSet[i]);
			}
		}
	}
	
	public void printInformativeKmers(int percCutoff){
		for(int i=0; i<MapMatrix.length; i++){
			int[][] kmerMap = MapMatrix[i];
			boolean pass = false;
			for(int j=0; j<kmerMap[0].length; j++){
				int colSum =0;
				for(int k=0; k<kmerMap.length; k++){
					colSum = colSum + kmerMap[j][k];
				}
				if((int)((colSum*100)/kmerMap.length) > percCutoff){
					pass = true;
				}
			}
			if(pass){
				System.out.println(Utilities.int2seq(i, k));
			}
			
		}
	}
	
	
	public static void main(String[] args){
		
		
	}

}
