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
	private double pbMax=2;
	
	public ChipSeq5PrimeProfiler(BinningParameters ps, SeqExpander exp, char strand, double pbMax) {
		params = ps;
		expanders = new ArrayList<SeqExpander>(); 
		expanders.add(exp);
		this.strand = strand;
	}
	public ChipSeq5PrimeProfiler(BinningParameters ps, List<SeqExpander> exps, char strand, double pbMax) {
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
			double[] exparray = new double[params.getNumBins()];
			for(int i = 0; i < exparray.length; i++) { exparray[i] = 0; }
			
			while(hits.hasNext()){
				SeqHit hit = hits.next();
				if (hit.getStrand()==this.strand){  //only count one strand
					if ((start<=hit.getFivePrime() && this.strand=='+')
							||(end>hit.getFivePrime() && this.strand=='-')){
						int hit5Prime = hit.getFivePrime()-start;
						exparray[params.findBin(hit5Prime)]+=hit.getWeight();
					}
				}				
			}
			for(int i = 0; i < array.length; i++) { 
				if(exparray[i]<=pbMax)
					array[i] += exparray[i];
				else
					array[i] += pbMax;
			}
		}		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	
	
	//No cleanup
	public void cleanup(){
		for(SeqExpander e : expanders)
			e.close();
	}
}
