/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.metaprofile;

import java.util.*;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;

public class ChipSeqProfiler implements PointProfiler<Point,PointProfile> {
	
	private Genome genome;
	private ExperimentManager manager=null;
	private BinningParameters params;
	private int extension; 
	private boolean useFivePrime=false;
	private char readStrand ='/';
	
	public ChipSeqProfiler(Genome gen, BinningParameters ps, ExperimentManager man, int ext, char strand) {
		genome = gen;
		manager = man;
		params = ps;
		extension=ext;
		if(extension==-1){
			useFivePrime = true;
			extension=0;
		}
		readStrand = strand;
	}
	
	public BinningParameters getBinningParameters() {
		return params;
	}
	
	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;
		
		
		int start = Math.max(1, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom()));
		
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		int ext = 200;
		Region extQuery = new Region(a.getGenome(), a.getChrom(), start-ext>0 ? start-ext : 1, end+ext < a.getGenome().getChromLength(a.getChrom()) ? end+ext : a.getGenome().getChromLength(a.getChrom()) );
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		for(ControlledExperiment expt : manager.getReplicates()){
			List<StrandedBaseCount> sbcs = expt.getSignal().getBases(extQuery);
			for(StrandedBaseCount sbc : sbcs){
				SeqHit hit = new SeqHit(genome, a.getChrom(), sbc);
				if(extension>0)
					hit = hit.extendHit(extension);
				
				if(hit.overlaps(query) && (readStrand=='.' || hit.getStrand()==readStrand)){
					int startOffset = Math.max(0, hit.getStart()-start);
					int endOffset = Math.max(0, Math.min(end, hit.getEnd()-start));
				
					if(!strand) { 
						int tmpEnd = window-startOffset;
						int tmpStart = window-endOffset;
						startOffset = tmpStart;
						endOffset = tmpEnd;
					}
					
					int startbin = params.findBin(startOffset);
					int endbin = params.findBin(endOffset);
					
					addToArray(startbin, endbin, array, 1.0);
				}	
			}
		}
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}

	private void addToArray(int i, int j, double[] array, double value) { 
		for(int k = i; k <= j; k++) { 
			array[k] += value;
		}
	}

	public void cleanup() {
	}
	
}
