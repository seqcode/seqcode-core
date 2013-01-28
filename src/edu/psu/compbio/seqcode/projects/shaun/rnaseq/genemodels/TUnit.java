package edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
/*
 * TUnit is a transcribed unit such as an exon or intron
 */
public class TUnit implements Comparable<TUnit> {
	protected Region coords;
	protected char strand='.';
	protected double hits=0;
	
	public TUnit(Region r, char str){
		coords = r;
		strand=str;
	}
	
	//Accessors
	public Region getCoords(){return coords;}
	public char getStrand(){return strand;}
	public double getHits(){return hits;}
	public void setCoords(Region r){coords=r;}
	
	public void addHit(){addHits(1.0);}
	public void addHits(double x){hits+=x;}

	public int compareTo(TUnit u) {
		return(coords.compareTo(u.coords));
	}
}
