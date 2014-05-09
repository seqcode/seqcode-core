package edu.psu.compbio.seqcode.gse.viz.metaprofile;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.projects.shaun.ConsensusSequence;
import edu.psu.compbio.seqcode.projects.shaun.ConsensusSequenceScoreProfile;
import edu.psu.compbio.seqcode.projects.shaun.ConsensusSequenceScorer;

public class ConsensusProfiler  implements PointProfiler<Point, Profile>{

	private ConsensusSequence consensus;
	private ConsensusSequenceScorer scorer;
	private SequenceGenerator seqgen;
	private Genome gen;
	private BinningParameters params=null;
	private double mismatchThreshold=0;
	private char searchStrand = '.'; //W=watson, C=crick. 
	
	public ConsensusProfiler(BinningParameters bp, Genome g, ConsensusSequence cons, double maxMismatch, boolean useCache, String seqPath, char watsoncrick){
		mismatchThreshold=maxMismatch;
		gen=g;
		params=bp; 
		consensus=cons;
		scorer = new ConsensusSequenceScorer(consensus);
		seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.setGenomePath(seqPath);
		}
	}

	public BinningParameters getBinningParameters() {
		return params;
	}

	public Profile execute(Point a) {
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		int start = Math.max(1, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom()));
		Region query = new Region(gen, a.getChrom(), start, end);
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;
		
		String seq = seqgen.execute(query);
		
		char strandChar = '.';
        if(searchStrand!='.'){
	        if(a instanceof StrandedPoint){
	        	char rstrand = ((StrandedPoint)a).getStrand(); 
	        	if((searchStrand=='W' && rstrand=='+') || (searchStrand=='C' && rstrand=='-'))
	        		strandChar='+';
	        	else
	        		strandChar='-';
	        }else
	        	strandChar = searchStrand=='W' ? '+' : '-';
        }
		ConsensusSequenceScoreProfile profiler = scorer.execute(seq, strandChar);
		for(int i=query.getStart(); i<query.getEnd(); i+=params.getBinSize()){
			boolean hasMatch=false;
			int matchOffset=-1;
			for(int j=i; j<i+params.getBinSize() && j<query.getEnd(); j++){
				int offset = j-query.getStart();
				if(profiler.getLowestMismatch(offset)<=mismatchThreshold){
					hasMatch=true;
					matchOffset = offset;
				}
			}
			if(hasMatch){
				int bin = params.findBin(matchOffset);
				if(!strand) { 
					bin = (params.getNumBins()-bin)-1;
				}
				maxToArray(bin, bin, array, 1);
			}
		}
		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}

	private void addToArray(int i, int j, double[] array, double value) { 
		for(int k = i; k <= j; k++) { 
			array[k] += value;
		}
	}
	private void maxToArray(int i, int j, double[] array, double value) { 
		for(int k = i; k <= j; k++) { 
			array[k] = Math.max(array[k],value);
		}
	}
	public void cleanup() {}
}