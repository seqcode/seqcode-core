package edu.psu.compbio.seqcode.gse.viz.metaprofile;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHitPair;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.projects.gps.ReadHit;

public class PairedSeqMidpointProfiler  implements PointProfiler<Point,PointProfile> {
	private BinningParameters params;
	private List<SeqExpander> expanders=null;
	private double perBaseMax=1000;
	
	public PairedSeqMidpointProfiler(BinningParameters ps, SeqExpander exp, double pbMax, char strand) {
		params = ps;
		expanders = new ArrayList<SeqExpander>(); 
		expanders.add(exp);
		perBaseMax = pbMax;
	}
	public PairedSeqMidpointProfiler(BinningParameters ps, List<SeqExpander> exps, double pbMax, char strand) {
		params = ps;
		expanders = exps;
		perBaseMax=pbMax;
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
		
		if(expanders!=null){
			for(SeqExpander expander : expanders){
				Iterator<SeqHitPair> pairs = expander.getPairs(extQuery);
				HashMap<Point, Double> readFilter = new HashMap<Point, Double>();
				
				while(pairs.hasNext()) {
					SeqHitPair pair = pairs.next();
					Point midpoint = pair.getMidpoint();
					
					if(query.contains(midpoint)){
						if(!readFilter.containsKey(midpoint))
							readFilter.put(midpoint, pair.getPairWeight());
						else
							readFilter.put(midpoint, readFilter.get(midpoint)+pair.getPairWeight());
						
						if(readFilter.get(midpoint)<=perBaseMax){
							int offset = Math.max(0, midpoint.getLocation()-start);
							
							if(!strand) { 
								int tmp = window-offset;
								offset = tmp;
							}
							
							int startbin = params.findBin(offset);
							
							array[startbin] += 1.0;
						}
					}
				}
			}
		}
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	public void cleanup(){
		if(expanders!=null){
			for(SeqExpander e : expanders)
				e.close();
		}
	}

}
