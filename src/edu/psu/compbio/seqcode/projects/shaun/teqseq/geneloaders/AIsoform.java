package edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders;

import java.util.ArrayList;
import java.util.Collections;


//ADDED BY YP
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ExptCollection;

public class AIsoform extends ARegion{

	protected ArrayList<ARegion> exons = new ArrayList<ARegion>();
	protected ArrayList<ARegion> CDSs = new ArrayList<ARegion>();
	protected ArrayList<ARegion> startCodons = new ArrayList<ARegion>();
	protected ArrayList<ARegion> stopCodons = new ArrayList<ARegion>();
	protected ARegion region5prime=null, region3prime=null;	//NOTE: These will stay null if it is Cufflinks type unless there is a strand.
	protected double[] coverage=null;
	protected StrandedPoint TSS = null;
	protected A3UTR utr = null;
		
	public AIsoform(Region c, char str, String name, String id, String transType, ARegionType regionType) {
		super(c, str, name, id, transType, regionType);
	}

	public ArrayList<ARegion> getExons(){return(exons);}
	public ArrayList<ARegion> getCDSs(){return(CDSs);}
	public ArrayList<ARegion> getStartCodons(){return(startCodons);}
	public ArrayList<ARegion> getStopCodons(){return(stopCodons);}
	public ARegion getRegion3Prime(){return region3prime;}
	public ARegion getRegion5Prime(){return region5prime;}
	public StrandedPoint getTSS(){return TSS;}
	public double getCoverage(int r){return (coverage==null || r>coverage.length || r<0)? 0.0 : coverage[r];}
	public void addCoverage(int r, double c){coverage[r]+=c;}
	
	public void initializeCoverageCounter(int reps){
		coverage = new double[reps];
		for(int r=0; r<coverage.length; r++){coverage[r]=0;}
	}
	
	public void addPart(ARegion r){
		switch (r.getRegionType()){
			case exon: exons.add(r); this.extendCoords(r.coords); break;
			case CDS: CDSs.add(r); this.extendCoords(r.coords); break;
			case start_codon: startCodons.add(r); this.extendCoords(r.coords); break;
			case stop_codon: stopCodons.add(r); this.extendCoords(r.coords); break;
		}
		
		if(this.getStrand()=='+' && (TSS==null || r.getCoords().getStart()<TSS.getLocation())){
			TSS = new StrandedPoint(coords.getGenome(), coords.getChrom(), r.getCoords().getStart(), strand);
		}else if(this.getStrand()=='-' && (TSS==null || r.getCoords().getEnd()>TSS.getLocation())){
			TSS = new StrandedPoint(coords.getGenome(), coords.getChrom(), r.getCoords().getEnd(), strand);
		}
		
		if(this.getStrand()=='+' && r.getRegionType()==ARegionType.exon && (region3prime==null || r.compareTo(region3prime)>0)){ region3prime=r;}
		if(this.getStrand()=='-' && r.getRegionType()==ARegionType.exon && (region3prime==null || r.compareTo(region3prime)<0)){ region3prime=r;}
		if(this.getStrand()=='+' && r.getRegionType()==ARegionType.exon && (region5prime==null || r.compareTo(region5prime)<0)){ region5prime=r;}
		if(this.getStrand()=='-' && r.getRegionType()==ARegionType.exon && (region5prime==null || r.compareTo(region5prime)>0)){ region5prime=r;}
	}
	
	public void sortRegions(){
		Collections.sort(exons);
		Collections.sort(CDSs);
		Collections.sort(startCodons);
		Collections.sort(stopCodons);
	}
	
	public String toString(){
		if (region5prime != null && region3prime != null){
			//System.out.println("Cufflinks");
			return(String.format(getName()+"\t"+getID()+"\t"+getTransType()+"\t"+getCoords().getLocationString()+":"+getStrand()+"\t"+exons.size()+" exons"+"\t"+CDSs.size()+" CDSs"+"\t"+startCodons.size()+" startCodons"+"\t"+stopCodons.size()+" stopCodons"));
		}
		return(String.format(getName()+"\t"+getID()+"\t"+getTransType()+"\t"+getCoords().getLocationString()+":"+getStrand()+"\t"+exons.size()+" exons"+"\t"+CDSs.size()+" CDSs"+"\t"+startCodons.size()+" startCodons"+"\t"+stopCodons.size()+" stopCodons"+"\n\t\t5':\t"+region5prime.getCoords()+"\n\t\t3':\t"+region3prime.getCoords()));
	}
	
	public void addUTR(A3UTR r){
		utr = r;
	}
	
	public A3UTR getUTR(){return utr;}
	////////////////////////////////////////
}
