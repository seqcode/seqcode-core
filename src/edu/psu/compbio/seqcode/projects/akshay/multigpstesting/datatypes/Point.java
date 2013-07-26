package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

import java.util.HashMap;
import java.util.Map;

public class Point {
	public String chr;
	public int coord;
	public static final int scanwidth = 500;
	public static final int infdistance = 201;
	public Tuple overlaps;
	public Point(String currentpoint){
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
	
	public Point(String chr, int coord, Tuple overlaps){
		this.chr = chr;
		this.coord = coord;
		this.overlaps = overlaps;
	}
	
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		
		Point pt = (Point) obj;
		return this.chr.equals(pt.chr) && this.coord == pt.coord;
	}
	
	@Override
	public int hashCode(){
		int result = 17;
		int code = (int) this.coord;
		code += (int) (this.chr == null ? 0 : this.chr.hashCode());
		result = 37*result + code;
		return result;
	}
	
	
	
	public int getDistancetoPoint(Point newPoint){
		if(newPoint.chr.matches(this.chr)){
			return Math.abs(newPoint.coord-this.coord);
		}
		else{
			return infdistance;
		}
	}
	
	public static void main(String [] args){
		Point pt1 = new Point("2:345");
		Point pt2 = new Point("2:344");
		System.out.println(pt1.equals(pt2));
		
		
		Map<Point,String> listOfPoints = new HashMap<Point,String>();
		listOfPoints.put(pt1, "rrr");
		System.out.println(listOfPoints.containsKey(pt2));
		
	}
	
	
}
