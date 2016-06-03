package org.seqcode.projects.shaun.teqseq.geneloaders;

import java.util.ArrayList;
import java.util.Collections;

import org.seqcode.genome.location.Region;
import org.seqcode.projects.shaun.teqseq.core.ExptCollection;


public class AGene extends ARegion{

	protected ArrayList<AIsoform> isoforms = new ArrayList<AIsoform>();
	protected ArrayList<A3UTR> utrs = new ArrayList<A3UTR>(); 
	protected ArrayList<A3UTR> novoUTR = new ArrayList<A3UTR>(); 
	protected double[] coverage=null;
	protected int coord3prime = -1; //The 3' most end of the CDS region for the gene
	
	public AGene(Region c, char str, String name, String id, String transType, ARegionType regionType){
		super(c, str, name, id, transType, regionType);
	}
	
	//Accessors
	public ArrayList<AIsoform> getIsoforms(){return(isoforms);}
	public double getCoverage(int r){return (coverage==null || r>coverage.length || r<0)? 0.0 : coverage[r];}
	public void addCoverage(int r, double c){coverage[r]+=c;}
	
	public void initializeCoverageCounter(int reps){
		coverage = new double[reps];
		for(int r=0; r<coverage.length; r++){coverage[r]=0;}
	}
	
	public void addIsoform(AIsoform i){
		isoforms.add(i);
		this.extendCoords(i.coords);
	}
	
	public void sortIsoforms(){
		for(AIsoform a : isoforms)
			a.sortRegions();
		Collections.sort(isoforms);
	}
	public String toString(){
		return(String.format(getName()+"\t"+getID()+"\t"+getTransType()+"\t"+getCoords().getLocationString()+":"+getStrand()+"\t"+isoforms.size()+" isoforms"));
	}
	
	public void set3coord(int coord){
		coord3prime = coord;
	}
	public int get3coord(){return coord3prime;}
	
	public void addUTR(A3UTR r){
		utrs.add(r);
	}
	public void addADD(A3UTR reg){utrs.add(reg);}
	public ArrayList<A3UTR> getUTR(){
		return utrs;
	}
	public void addNovoUTR(A3UTR r){
		novoUTR.add(r);
	}
	public ArrayList<A3UTR> getNovoUTR(){
		return novoUTR;
	}
	public void setNovoUTR(ArrayList<A3UTR> utrs){
		novoUTR = utrs;
	}
	public void removeNovoUTR(int index){
		novoUTR.remove(index);
	}
	public int get3End(){
		int endPos = -1;
		int endNeg = Integer.MAX_VALUE;
		int end = 0;
		
		if (getStrand()=='+'){
			for (ARegion annot : utrs){
				if (annot.getCoords().getEnd() >= endPos){endPos = annot.getCoords().getEnd();}
			}
			for (ARegion novo : novoUTR){
				if (novo.getCoords().getEnd() >= endPos){endPos = novo.getCoords().getEnd();}
			}
			end = endPos;
		}
		else if (getStrand()=='-'){
			for (ARegion annot : utrs){
				if (annot.getCoords().getStart() <= endNeg){endNeg = annot.getCoords().getStart();}
			}
			for (ARegion novo : novoUTR){
				if (novo.getCoords().getStart() <= endNeg){endNeg = novo.getCoords().getStart();}
			}
			end = endNeg;
		}
		return end;
	}
	public String toStringGFF3(){
		String result = new String();
		if (getStrand()=='+'){
			result = String.format("chr"+getCoords().getChrom()+"\t"+getTransType()+"\t"+"gene"+"\t"+coord3prime+"\t"+get3End()+"\t"+"."+"\t"+getStrand()+"\t"+"."+"\t"+"ID="+getID()+";"+"Name="+getName()+"\n");
		}
		else if (getStrand()=='-'){
			result = String.format("chr"+getCoords().getChrom()+"\t"+getTransType()+"\t"+"gene"+"\t"+get3End()+"\t"+coord3prime+"\t"+"."+"\t"+getStrand()+"\t"+"."+"\t"+"ID="+getID()+";"+"Name="+getName()+"\n");
		}
		
		for (A3UTR utr : getUTR()){
			if (utr != null){
				result = result + String.format("chr"+getCoords().getChrom()+"\t"+getTransType()+"\t"+"mRNA"+"\t"+utr.getCoords().getStart()+"\t"+utr.getCoords().getEnd()+"\t"+"."+"\t"+getStrand()+"\t"+"."+"\t"+"ID="+utr.getID()+";Name="+getName()+";Parent="+getID()+"\n");
				int index = 1; //since the reg's ID would have exon_# in it, must add it in by hand
				for (ARegion reg : utr.getRegions()){
					result = result + String.format("chr"+getCoords().getChrom()+"\t"+getTransType()+"\t"+"exon"+"\t"+reg.getCoords().getStart()+"\t"+reg.getCoords().getEnd()+"\t"+"."+"\t"+getStrand()+"\t"+"."+"\t"+"ID="+utr.getID()+"."+index+";Name="+getName()+";Parent="+utr.getID()+"\n");
					index++;
				}
			}
		}
		
		for (A3UTR novo : novoUTR){
			result = result + String.format("chr"+getCoords().getChrom()+"\t"+getTransType()+"\t"+"mRNA"+"\t"+novo.getCoords().getStart()+"\t"+novo.getCoords().getEnd()+"\t"+"."+"\t"+getStrand()+"\t"+"."+"\t"+"ID="+novo.getID()+";Name="+getName()+";Parent="+getID()+"\n");
			int index = 1; //since the reg's ID would have exon_# in it, must add it in by hand
			for (ARegion exon : novo.getRegions()){
				result = result + String.format("chr"+getCoords().getChrom()+"\t"+getTransType()+"\t"+"exon"+"\t"+exon.getCoords().getStart()+"\t"+exon.getCoords().getEnd()+"\t"+"."+"\t"+getStrand()+"\t"+"."+"\t"+"ID="+novo.getID()+"."+index+";Name="+getName()+";Parent="+novo.getID()+"\n");
				index++;
			}
		}
		return result;
	}
}
