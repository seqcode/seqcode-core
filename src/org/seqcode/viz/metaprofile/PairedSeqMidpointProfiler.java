package org.seqcode.viz.metaprofile;

import java.util.List;

import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;


public class PairedSeqMidpointProfiler  implements PointProfiler<Point,PointProfile> {
	private Genome genome;
	private BinningParameters params;
	private ExperimentManager manager=null;
	private double perBaseMax=1000;
	
	public PairedSeqMidpointProfiler(Genome gen, BinningParameters ps, ExperimentManager man) {
		genome = gen;
		manager = man;
		params = ps;
		/*******NOTE: ExperimentManager needs to be able to handle paired reads before any of this will work****/
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
			List<StrandedPair> sps = expt.getSignal().getPairs(extQuery);
			for(StrandedPair pair : sps){
				Point midpoint = pair.getMidpoint();
				if(pair!=null){
					if(query.contains(midpoint)){
						if (start<=midpoint.getLocation() && end>=midpoint.getLocation()){
								int hitLoc = midpoint.getLocation()-start;
								if(!strand)
									hitLoc = end-midpoint.getLocation();
								array[params.findBin(hitLoc)]+=pair.getWeight();
						}
					}
				}
			}
		}
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	public void cleanup(){
	}

}
