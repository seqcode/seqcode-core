/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.metagenes;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.*;

public class ChipSeq5PrimeProfiler implements PointProfiler<Point,PointProfile> {
	
	private BinningParameters params;
	private List<SeqExpander> expanders;
	private char strand;
	
	public ChipSeq5PrimeProfiler(BinningParameters ps, SeqExpander exp, char strand) {
		params = ps;
		expanders = new ArrayList<SeqExpander>(); 
		expanders.add(exp);
		this.strand = strand;
	}
	public ChipSeq5PrimeProfiler(BinningParameters ps, List<SeqExpander> exps, char strand) {
		params = ps;
		expanders = exps;
		this.strand = strand;
	}

	public BinningParameters getBinningParameters() {
		return params;
	}

	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		for(SeqExpander expander : expanders){
			Iterator<SeqHit> hits = expander.execute(query);
			List <SeqHit> hitList =  filterDuplicateHits(hits);

			for(SeqHit hit : hitList) {
				if (hit.getStrand()==this.strand){  //only count one strand
					if ((start<=hit.getFivePrime() && this.strand=='+')
							||(end>hit.getFivePrime() && this.strand=='-')){
						int hit5Prime = hit.getFivePrime()-start;
						array[params.findBin(hit5Prime)]++;
					}
				}				
			}
		}		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	/*
	 * filter out duplicate reads (potential tower, needles, but could be real reads)
	 * assuming the reads are sorted
	 */
	public List<SeqHit> filterDuplicateHits(Iterator<SeqHit> hits){
		SeqHit currentHit = hits.next();
		int count=1;
		List<SeqHit> filteredReads = new ArrayList<SeqHit>();
		while(hits.hasNext()) {
			SeqHit hit = hits.next();
			// if read from a new position
			if (!(currentHit.getStart()==hit.getStart())){
				currentHit = hit;
				count=0;
				filteredReads.add(hit);
			}
			else {// if  duplicate
				count++;
				if (count<=3)
					filteredReads.add(hit);
			}
		}
		return filteredReads;
	}
	
	//No cleanup
	public void cleanup(){
		for(SeqExpander e : expanders)
			e.close();
	}
}
