package org.seqcode.data.seqdata;

import org.seqcode.genome.location.Point;

public class SeqHitPair {

	private SeqHit left, right;
	private double pairWeight = 1.0;
	private int pairCode;

	public SeqHitPair(SeqHit l, SeqHit r, double pweight, int pcode) {
		left = l;
		right = r;
		pairWeight = pweight;
		pairCode = pcode;
	}
	
	public SeqHit getLeft(){return left;}
	public SeqHit getRight(){return right;}
	
	public double getPairWeight() {	return pairWeight; }
	public int getCode(){ return pairCode; }
		
	public void setPairWeight(double weight) { this.pairWeight = weight;}
	
	public Point getMidpoint(){
		if(left.getChrom() != right.getChrom())
			return null;
		else{
			return new Point(left.getGenome(), left.getChrom(), (left.getFivePrime()+right.getFivePrime())/2);
		}
	}
		
}
