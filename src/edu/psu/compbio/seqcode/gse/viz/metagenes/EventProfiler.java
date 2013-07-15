package edu.psu.compbio.seqcode.gse.viz.metagenes;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class EventProfiler implements PointProfiler<Point, Profile>{

	private Genome gen;
	private BinningParameters params=null;
	private List<Point> events;
	
	public EventProfiler(BinningParameters bp, Genome g, ArrayList<Point> e){
		gen=g;
		params=bp; 
		events = e;
		System.out.println(events.size());
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
		//candidates
		ArrayList<Point> candidates = new ArrayList<Point>();
		for(Point p : events){
			if(query.contains(p))
				candidates.add(p);
		}
		
		for(int i=query.getStart(); i<query.getEnd(); i+=params.getBinSize()){
			Region curr = new Region(query.getGenome(), query.getChrom(), i, i+params.getBinSize());
			for(Point p : candidates){
				if(curr.contains(p)){
					int offset = p.getLocation()-query.getStart();
					if(!strand) { 
						int tmp = window-offset;
						offset = tmp;
					}
					int bin = params.findBin(offset);
					addToArray(bin, bin, array, 1);
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
	public void cleanup() {}
}
