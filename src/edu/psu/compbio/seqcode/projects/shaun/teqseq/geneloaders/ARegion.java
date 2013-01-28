package edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
/*
 * ARegion = AnnotatedRegion. 
 * 
 */
public class ARegion implements Comparable<ARegion>{
	
	protected String name;
	protected String ID;
	protected String transType;
	protected ARegionType regionType;
	protected Region coords; //Outline of the annotated region
	protected char strand = '.';
	
	public ARegion(Region c, char str, String name, String id, String transType, ARegionType regionType){
		coords = c;
		strand = str;
		this.name = name;
		this.ID = id;
		this.transType = transType;
		this.regionType = regionType;
	}
	
	//Accessors
	public Region getCoords(){return coords;}
	public char getStrand(){return strand;}
	public String getName(){return name;}
	public String getID(){return ID;}
	public String getTransType(){return transType;}
	public ARegionType getRegionType(){return regionType;}
	
	//Extend the current coordinates (extends only, never contracts)
	public void extendCoords(Region newCoords){
		if(!coords.contains(newCoords)){
			//Two ways to interpret this. 
			if(newCoords.contains(coords))//If the new coordinates fully contain the old, this is an expanded coordinate set
				coords = newCoords;
			else{ //Otherwise, this is an exon that is being added to the gene. expand the coordinates in that direction
				int dstart=0, dend=0;
				if(newCoords.getStart()<coords.getStart())
					dstart = coords.getStart()-newCoords.getStart();
				if(newCoords.getEnd()>coords.getEnd())
					dend = newCoords.getEnd() - coords.getEnd();
				extendCoords(dstart, dend);
			}
		}
	}
	public void extendCoords(int startExt, int endExt){
		coords = coords.expand(startExt, endExt);
	}
	
	public static ARegionType translateToARegionType(String x){
		if(x.equals("exon"))
			return(ARegionType.exon);
		else if(x.equals("CDS"))
			return(ARegionType.CDS);
		else if(x.equals("start_codon"))
			return(ARegionType.start_codon);
		else if(x.equals("stop_codon"))
			return(ARegionType.stop_codon);
		else
			return(null);
	}
	
	public String toString(){
		return(String.format(getName()+"\t"+getID()+"\t"+getTransType()+"\t"+getCoords().getLocationString()+":"+getStrand()));
	}

	public int compareTo(ARegion o) {
		return(this.coords.compareTo(o.coords));
	}
}
