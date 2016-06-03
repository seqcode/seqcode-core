package org.seqcode.gse.viz.metaprofile;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.utils.sequence.SequenceUtils;
import org.seqcode.projects.shaun.ConsensusSequence;
import org.seqcode.projects.shaun.ConsensusSequenceScoreProfile;
import org.seqcode.projects.shaun.ConsensusSequenceScorer;

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
		Region query = a.expand(window/2);
		
		char rstrand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() : '.';
		String seq = seqgen.execute(query);
		if(rstrand=='-')
			seq = SequenceUtils.reverseComplement(seq);
		
		ConsensusSequenceScoreProfile profiler = scorer.execute(seq, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
		for(int i=0; i<seq.length() && i<window; i+=params.getBinSize()){
			if(profiler.getLowestMismatch(i)<=mismatchThreshold){
				int bin = params.findBin(i);
				addToArray(bin, bin, array, 1);
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