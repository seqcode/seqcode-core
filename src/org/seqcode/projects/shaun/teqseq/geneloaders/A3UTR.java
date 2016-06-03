package org.seqcode.projects.shaun.teqseq.geneloaders;

import java.util.ArrayList;
import java.util.Collections;

import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.projects.shaun.teqseq.geneloaders.ARegion;
import org.seqcode.projects.shaun.teqseq.geneloaders.ARegionType;


public class A3UTR extends ARegion {
	
	public ArrayList<ARegion> regions = new ArrayList<ARegion>();

	public A3UTR(Region c, char str, String name, String id, String transType, ARegionType regionType) {
		super(c, str, name, id, transType, regionType);
	}
	
	public ArrayList<ARegion> getRegions(){return regions;}
	public void inputRegions(ArrayList<ARegion> rs){regions = rs;}
	
	public void filterRegions(ArrayList<ARegion> rs, int endCoord){ //endCoord is the 3'most CDS end coordinate
		ArrayList<ARegion> set = new ArrayList<ARegion>();
		for (ARegion reg : rs){
			int start = reg.getCoords().getStart();
			int end = reg.getCoords().getEnd();
			
			if (getStrand()=='+'){
				if (start > endCoord){set.add(reg);}
				else if (endCoord>=start && endCoord<end){
					Region r = new Region(getCoords().getGenome(), getCoords().getChrom(), endCoord, end);
					ARegion areg = new ARegion(r, getStrand(), getName(), getID(), getTransType(), ARegionType.exon);
					set.add(areg);
				}
			}
			else if (getStrand()=='-'){
				if (end < endCoord){set.add(reg);}
				else if (endCoord>start && endCoord<=end){
					Region r = new Region(getCoords().getGenome(), getCoords().getChrom(), start, endCoord);
					ARegion areg = new ARegion(r, getStrand(), getName(), getID(), getTransType(), ARegionType.exon);
					set.add(areg);
				}
			}
		}
		Collections.sort(set);
		if (getStrand()=='+' && set.get(0).getCoords().getStart()!=endCoord){
			int end = set.get(0).getCoords().getEnd();
			Region r = new Region(getCoords().getGenome(), getCoords().getChrom(), endCoord, end);
			set.set(0, new ARegion(r, getStrand(), getName(), getID(), getTransType(), ARegionType.exon));
		}
		else if (getStrand()=='-' && set.get(set.size()-1).getCoords().getEnd()!=endCoord){
			int start = set.get(set.size()-1).getCoords().getStart();
			Region r = new Region(getCoords().getGenome(), getCoords().getChrom(), start, endCoord);
			set.set(set.size()-1, new ARegion(r, getStrand(), getName(), getID(), getTransType(), ARegionType.exon));
		}
		regions = set;
	}
	
	//add regions
	public void addRegion(ARegion r){
		if (r != null){
			regions.add(r);
			this.extendCoords(r.getCoords());
		}
	}
	
	//toString
	public String toString(){
		return(String.format(getName()+"\t"+getID()+"\t"+getTransType()+"\t"+getCoords().getLocationString()+":"+getStrand()+"\t"+regions.size()+" exons spanned"));
	}
}
