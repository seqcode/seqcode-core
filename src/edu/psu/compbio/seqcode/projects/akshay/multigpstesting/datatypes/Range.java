package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

public class Range {
	private int lowerBound;
	private int upperBound;
	private String chr;
		
	public Range(String chr, int lowerBound, int upperBound){
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.chr = chr;
	}
		
	boolean includes(Point givenPoint){
		if(givenPoint.chr.matches(this.chr)){
			return givenPoint.coord >= this.lowerBound && givenPoint.coord < this.upperBound;
		}
		else{
			return false;
		}
			
	}
		
	boolean includes(String givenPoint){
		String[] pieces = givenPoint.split(":");
		if (pieces[0].matches(chr)){
			return Integer.parseInt(pieces[1]) >= lowerBound && Integer.parseInt(pieces[1]) < upperBound;
		}
		else{
			return false;
		}
	}
		
	String getChr(){
		return chr;
	}

}


