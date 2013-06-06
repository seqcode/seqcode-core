/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.metagenes;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.*;

public class SimpleChipSeqProfilerBAM implements PointProfiler<Point,PointProfile> {
	
	private BinningParameters params;
	private DeepSeqExpt expt;
	private int extension; 
	private double perBaseMax=100;
	private boolean useFivePrime = false;
	private char readStrand ='/';
	
	public SimpleChipSeqProfilerBAM(BinningParameters ps, DeepSeqExpt exp) { 
		this(ps, exp, 175, 100, '/');
	}
	public SimpleChipSeqProfilerBAM(BinningParameters ps, DeepSeqExpt exp, int ext, double pbMax, char strand) {
		params = ps;
		expt = exp;
		extension=ext;
		if(extension==-1)
			useFivePrime=true;
		perBaseMax = pbMax;
		readStrand = strand;
	}
	

	public BinningParameters getBinningParameters() {
		return params;
	}
	public void setUseFivePrime(boolean ufp){useFivePrime = ufp;}

	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;
		
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		Region extQuery = new Region(a.getGenome(), a.getChrom(), start-extension>0 ? start-extension : 1, end+extension < a.getGenome().getChromLength(a.getChrom()) ? end+extension : a.getGenome().getChromLength(a.getChrom()) );
		
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		List<ReadHit> hits = expt.loadHits(extQuery);
		HashMap<Region, Double> readFilter = new HashMap<Region, Double>();
		
		for(ReadHit hit : hits){
			if(hit.overlaps(query) && (readStrand=='/' || hit.getStrand()==readStrand)){
				if(!readFilter.containsKey(hit))
					readFilter.put(hit, hit.getWeight());
				else
					readFilter.put(hit, readFilter.get(hit)+hit.getWeight());
				
				if(readFilter.get(hit)<=perBaseMax){
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
	
	public void cleanup(){
		expt.closeLoaders();
	}
}
