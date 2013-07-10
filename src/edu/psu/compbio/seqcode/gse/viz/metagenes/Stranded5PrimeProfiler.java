/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.metagenes;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.*;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class Stranded5PrimeProfiler implements PointProfiler<Point,PointProfile> {
	
	private BinningParameters params;
	private List<SeqExpander> expanders=null;
	private DeepSeqExpt expt=null;
	private char strand;
	private double pbMax=2;
	
	public Stranded5PrimeProfiler(BinningParameters ps, SeqExpander exp, char strand, double pbMax) {
		params = ps;
		expanders = new ArrayList<SeqExpander>(); 
		expanders.add(exp);
		this.strand = strand;
	}
	public Stranded5PrimeProfiler(BinningParameters ps, List<SeqExpander> exps, char strand, double pbMax) {
		params = ps;
		expanders = exps;
		this.strand = strand;
	}
	public Stranded5PrimeProfiler(BinningParameters ps, DeepSeqExpt exp, char strand, double pbMax) {
		params = ps;
		expt = exp;
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
		char pointStrand = '+';
		if(a instanceof StrandedPoint)
			pointStrand = ((StrandedPoint)a).getStrand();
		boolean wantPosStrandReads = this.strand=='+';
		if(pointStrand == '-')
			wantPosStrandReads = !wantPosStrandReads;
		char wantedStrand = wantPosStrandReads?'+':'-';
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		double[] exparray = new double[params.getNumBins()];
		for(int i = 0; i < exparray.length; i++) { exparray[i] = 0; }
		
		if(expt!=null){
			Pair<ArrayList<Integer>, ArrayList<Float>> sbc = expt.loadStrandedBaseCounts(query, wantedStrand);
			for(int x=0; x<sbc.car().size(); x++){
				int pos = sbc.car().get(x);
				float weight = sbc.cdr().get(x);
				int hit5Prime = pos-start;
				exparray[params.findBin(hit5Prime)]+=weight;
			}
		}else if (expanders!=null){
			for(SeqExpander expander : expanders){
				Iterator<SeqHit> hits = expander.execute(query);
				while(hits.hasNext()){
					SeqHit hit = hits.next();
					if (hit.getStrand()==wantedStrand){  //only count one strand
						if (start<=hit.getFivePrime() && end>hit.getFivePrime()){
							int hit5Prime = hit.getFivePrime()-start;
							exparray[params.findBin(hit5Prime)]+=hit.getWeight();
						}
					}				
				}
			}
		}
		
		for(int i = 0; i < array.length; i++) { 
			if(exparray[i]>pbMax)
				exparray[i]=pbMax;
		}
		for(int i = 0; i < array.length; i++) {
			if(pointStrand == '+')
				array[i] += exparray[i];
			else
				array[i] += exparray[params.getNumBins()-i-1];
		}
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	
	//No cleanup
	public void cleanup(){
		if(expanders!=null){
			for(SeqExpander e : expanders)
				e.close();
		}
		if(expt!=null)
			expt.closeLoaders();
	}
}
