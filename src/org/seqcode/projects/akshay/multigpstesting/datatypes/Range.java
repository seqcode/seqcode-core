package org.seqcode.projects.akshay.multigpstesting.datatypes;

public class Range {
	public int lowerBound;
	public int upperBound;
	public String chr;
		
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
	
	@Override
	public boolean equals(Object obj){
		if (obj == this){
			return true;
		}
		
		Range ra = (Range) obj;
		return this.lowerBound == ra.lowerBound && this.upperBound == ra.upperBound && this.chr.equals(ra.chr);
	}
	
	@Override
	public int hashCode(){
		int result = 17;
		int code = (int) lowerBound;
		code += (int) upperBound;
		code += (int) (this.chr == null ? 0 :this.chr.hashCode());
		result = result*37 + code;
		return result;
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
	
	public String getRange(){
		return this.chr+":"+Integer.toString(this.lowerBound)+":"+Integer.toString(this.upperBound);
	}
	
	public static void main(String[] args){
		Range ra1 = new Range("4",450,550);
		Range ra2 = new Range("4",450,550);
		System.out.println(ra1.equals(ra2));
	}

}


