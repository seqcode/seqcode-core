package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

public class Point {
	public String chr;
	public int coord;
	public static final int scanwidth = 500;
	public static final int infdistance = 1000;
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
		
	public int getDistancetoPoint(Point newPoint){
		if(newPoint.chr.matches(this.chr)){
			return Math.abs(newPoint.coord-this.coord);
		}
		else{
			return infdistance;
		}
	}
}
