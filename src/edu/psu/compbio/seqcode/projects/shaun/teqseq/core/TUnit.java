package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;

public class TUnit implements Comparable<TUnit>{
	protected Region coords;
	protected char strand='.';
	protected int ID;
	protected double hitCount=0;
	protected double hitWeight=0;
	protected double[] hitStarts=null;
	protected double[] coverage=null;
	protected int numReplicates =0;
	protected List<String> annotations = new ArrayList<String>();
	protected List<AlignHit> boundarySpanningHits = null;
	protected List<PairedHit> boundarySpanningPairedHits = null;
	
	public TUnit(Genome g, String c, int s, int e, int i, char str){
		coords = new Region(g,c,s,e);	
		ID = i;
		strand = str;
	}
	public TUnit(Region reg, int i, char str){
		coords = reg;	
		ID = i;
		strand = str;
	}
	
	//Accessors
	public Region getCoords(){return(coords);}
	public int getID(){return ID;}
	public char getStrand(){return strand;}
	public int getRepCount(){return numReplicates;}
	public double getHitCount(){return hitCount;}
	public double getHitWeight(){return hitWeight;}
	public double getCoverage(int r){return (coverage==null || r>coverage.length || r<0)? 0.0 : coverage[r];}
	public double getHitStarts(int r){return (hitStarts==null || r>hitStarts.length || r<0)? 0.0 : hitStarts[r];}
	public void setCoords(Region r){coords = r;}
	public void addHit(double w){hitCount++; hitWeight+=w;}
	public void addCoverage(int r, double c){coverage[r]+=c;}
	public void addHitStart(int r, double h){hitStarts[r]+=h;}
	public List<AlignHit> getBoundaryHits(){return boundarySpanningHits;}
	public List<PairedHit> getBoundaryPairedHits(){return boundarySpanningPairedHits;}
	
	public void addBoundaryHit(AlignHit h){
		if(boundarySpanningHits == null)
			boundarySpanningHits = new ArrayList<AlignHit>();
		boundarySpanningHits.add(h);
	}
	public void addBoundaryPairedHit(PairedHit h){
		if(boundarySpanningPairedHits == null)
			boundarySpanningPairedHits = new ArrayList<PairedHit>();
		boundarySpanningPairedHits.add(h);
	}
	
	public void addAnnotation(String ID, String name, String postfix){
		if(annotations==null)
			annotations = new ArrayList<String>();
		annotations.add(ID+":"+name+":"+postfix);
	}
	
	public List<String> getAnnotations(){return annotations;}
	
	public String getAnnotationString(){
		if(annotations==null || annotations.size()==0){
			return "NONE";
		}
		String a ="";
		for(String x : annotations){
			a = a+x+";";
		}
		return(a);
	}
	
	public void initializeCounters(ExptCollection expts){initializeCounters(expts.getReplicateCount());}
	public void initializeCounters(int numRep){
		numReplicates = numRep;
		coverage = new double[numRep];
		hitStarts = new double[numRep];
		for(int r=0; r<coverage.length; r++){coverage[r]=0; hitStarts[r]=0;}
	}
	
	public boolean overlaps(AlignHit hit){
		return(coords.overlaps(hit));
	}
	
	public int compareTo(TUnit o) {
		return(this.coords.compareTo(o.coords));
	}
}
