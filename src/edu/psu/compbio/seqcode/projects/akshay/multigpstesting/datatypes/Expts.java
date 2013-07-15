package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

import java.util.*;


public class Expts {
	public Map<Point,String> isotaedPoints= new HashMap<Point,String>();
	public Map<Range,Integer> blacklist = new HashMap<Range,Integer>();
	
	public void mapPointstoRnages(Point currentPoint){
		if(blacklist.containsKey(currentPoint.overlaps.getFirst()){
			
		}
		
	}
	
	
}

class Point{
	public String chr;
	public int coord;
	public static final int scanwidth = 500; 
	public Tuple overlaps;
	Point(String currentpoint){
		String[] pieces = currentpoint.split(":");
		this.chr = pieces[0];
		this.coord = Integer.parseInt(pieces[1]);
		int quotient = this.coord/scanwidth;
		int reminder = this.coord%scanwidth;
		Range first = new Range(this.chr,scanwidth*quotient,scanwidth*quotient+scanwidth);
		Range last=null;
		if (reminder >= scanwidth/2){
			last = new Range(this.chr,scanwidth*quotient+scanwidth/2,scanwidth*quotient+scanwidth+scanwidth/2);
		}
		else{
			last = new Range(this.chr,scanwidth*quotient-scanwidth/2,scanwidth*quotient+scanwidth-scanwidth/2);
		}
		this.overlaps = new Tuple(first,last);
	}

}

class Tuple{
	public Range first;
	public Range last;
	Tuple(Range first, Range last){
		this.first = first;
		this.last = last;
	}
	Range getFirst(){
		return this.first;
	}
	Range getLast(){
		return this.last;
	}
}