package edu.psu.compbio.seqcode.gse.viz.metaprofile;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

/**
 * EventProfiler profiles the occurrence of point events around a set of reference points. 
 * 
 * @author mahony
 *
 */
public class EventProfiler implements PointProfiler<Point, Profile>{

	private Genome gen;
	private BinningParameters params=null;
	private Map<String, List<Point>> events; //events indexed by chromosome and then sorted
	
	public EventProfiler(BinningParameters bp, Genome g, List<Point> elist){
		gen=g;
		params=bp; 
		events = new HashMap<String, List<Point>>();
		
		for(Point e : elist){
			String currChr = e.getChrom();
			if(!events.containsKey(currChr)){
				events.put(currChr, new ArrayList<Point>());
			}
			events.get(currChr).add(e);
		}
		for(String c : events.keySet())
			Collections.sort(events.get(c));
	}

	public BinningParameters getBinningParameters() {
		return params;
	}

	/**
	 * execute: generate a PointProfile based on a new Point
	 */
	public Profile execute(Point a) {
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		//Initialize the query region
		int start = Math.max(1, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom()));
		Region query = new Region(gen, a.getChrom(), start, end);
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;

		//Find the candidate events using a binary search
		List<Point> candidates = new ArrayList<Point>();
		if(events.containsKey(a.getChrom())){
			List<Point> chrEvents = events.get(a.getChrom());
			int index = Collections.binarySearch(chrEvents, query.startPoint());
			for(int i=index; i<chrEvents.size(); i++){
				Point p = chrEvents.get(i);
				if(query.contains(p))
					candidates.add(p);
				else if(p.getLocation()>query.getEnd())
					break;
			}
		}
		
		//Make the new PointProfile
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
					array[bin]++;
				}
			}
		}
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}

	public void cleanup() {}
}
