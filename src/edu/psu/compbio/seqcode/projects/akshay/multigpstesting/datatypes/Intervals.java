package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

import java.util.*;


public class Intervals {
	public List<Range> ranges = new ArrayList<Range>();
	void add(int value, String chr, int lowerBound, int upperBound){
		this.ranges.add(new Range(chr,lowerBound,upperBound));
	}
	
	void add(Range range){
		ranges.add(range);
	}
	
	boolean inList(String givenPoint){
		boolean present = false;
		for(Range range: ranges){
			if (range.includes(givenPoint)){
				present = true;
			}
		}
		return present;
	}

	int getNoOfIntervals(){
		return ranges.size();
	}
}
