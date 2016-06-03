package org.seqcode.projects.shaun.rnaseq.genemodels;

import java.util.Iterator;
import java.util.TreeSet;

import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;

/*
 * SplicedTUnit defines a transcript. It is a TUnit and serves to hold other TUnits such as exons and introns
 */
public class SplicedTUnit extends TUnit{
	protected String name;
	protected String ID;
	protected String type;
	protected TreeSet<TUnit> exons = new TreeSet<TUnit>();
	protected TreeSet<TUnit> CDSs = new TreeSet<TUnit>();
	protected TreeSet<TUnit> introns = new TreeSet<TUnit>();
	protected Point TSS=null;
		
	public SplicedTUnit(String name, String ID, String type, Region r, char str) {
		super(r, str);
		this.name = name;
		this.ID = ID;
		this.type = type;
	}
	
	//Accessors
	public String getName(){return name;}
	public String getID(){return ID;}
	public String getType(){return type;}
	public int getNumExons(){return exons.size();}
	public Iterator<TUnit> getExonIterator(){return exons.iterator();}
	public Iterator<TUnit> getIntronIterator(){return introns.iterator();}
	public Iterator<TUnit> getCDSIterator(){return CDSs.iterator();}
	public Point getTSS(){return(TSS);}

	
	public void addExon(TUnit e){
		//Do we have this exon already? 
		boolean haveExon = false;
		for(TUnit u : exons){
			if(e.getCoords().matchesRegion(u.getCoords()))
				haveExon = true;
		}
		if(!haveExon){
			//Expand the region if exon lies outside
			if(!coords.contains(e.getCoords())) { 
				if(e.getCoords().getChrom().equals(coords.getChrom())){
					int upstream=0, downstream=0; 
					if(e.getCoords().getStart()<coords.getStart()){
						upstream = coords.getStart() - e.getCoords().getStart();
					}if(e.getCoords().getEnd()>coords.getEnd()){
						downstream = e.getCoords().getEnd()-coords.getEnd();
					}
					coords = coords.expand(upstream, downstream);
				}else{
					throw new IllegalArgumentException();
				}
			}
			//Add the exon
		    if(exons == null) { exons = new TreeSet<TUnit>(); }
	        exons.add(e);
	        consolidateExons();
	        findTSS();
	        findIntrons();
		}
	}
	
	public void addCDS(TUnit c){
		//Do we have this exon already? 
		boolean haveCDS = false;
		for(TUnit u : CDSs){
			if(c.getCoords().matchesRegion(u.getCoords()))
				haveCDS = true;
		}
		if(!haveCDS){
			//Expand the region if exon lies outside
			if(!coords.contains(c.getCoords())) { 
				if(c.getCoords().getChrom().equals(coords.getChrom())){
					int upstream=0, downstream=0; 
					if(c.getCoords().getStart()<coords.getStart()){
						upstream = coords.getStart() - c.getCoords().getStart();
					}if(c.getCoords().getEnd()>coords.getEnd()){
						downstream = c.getCoords().getEnd()-coords.getEnd();
					}
					coords = coords.expand(upstream, downstream);
				}else{
					throw new IllegalArgumentException();
				}
			}
			//Add the CDS
		    if(CDSs == null) { CDSs = new TreeSet<TUnit>(); }
	        CDSs.add(c);
		}
	}
	
	public void consolidateExons(){
		TreeSet<TUnit> conExons = new TreeSet<TUnit>();
		TreeSet<TUnit> conCDSs = new TreeSet<TUnit>();
		
		TUnit curr=null;
		Iterator<TUnit> exIter = exons.iterator();
		while(exIter.hasNext()){
		    TUnit e = exIter.next();
		    if(curr != null && e.getCoords().overlaps(curr.getCoords())){
		    	int upstream=0, downstream=0; 
				if(e.getCoords().getStart()<curr.getCoords().getStart()){
					upstream = curr.getCoords().getStart() - e.getCoords().getStart();
				}if(e.getCoords().getEnd()>curr.getCoords().getEnd()){
					downstream = e.getCoords().getEnd()-curr.getCoords().getEnd();
				}
		    	curr.setCoords(curr.getCoords().expand(upstream, downstream));
		    }else{
		    	curr=e;
		    	conExons.add(e);
		    }
		}
		
		curr=null;
		Iterator<TUnit> cdsIter = CDSs.iterator();
		while(cdsIter.hasNext()){
		    TUnit e = cdsIter.next();
		    if(curr != null && e.getCoords().overlaps(curr.getCoords())){
		    	int upstream=0, downstream=0; 
				if(e.getCoords().getStart()<curr.getCoords().getStart()){
					upstream = curr.getCoords().getStart() - e.getCoords().getStart();
				}if(e.getCoords().getEnd()>curr.getCoords().getEnd()){
					downstream = e.getCoords().getEnd()-curr.getCoords().getEnd();
				}
		    	curr.setCoords(curr.getCoords().expand(upstream, downstream));
		    }else{
		    	curr=e;
		    	conCDSs.add(e);
		    }
		}
		exons = conExons;
		CDSs = conCDSs;
		
		this.findIntrons();
	}
	
	public void findIntrons(){
		introns = new TreeSet<TUnit>();
		boolean first=true;
		TUnit lastE = null;
		Iterator<TUnit> exIter = exons.iterator();
		while(exIter.hasNext()){
		    TUnit e = exIter.next();
		    if(!first){
		    	Region currIntron = new Region(coords.getGenome(), coords.getChrom(), lastE.getCoords().getEnd()-1, e.getCoords().getStart()-1);
		    	introns.add(new TUnit(currIntron, e.getStrand()));
		    }
		    lastE = e;
		    first=false;
		}
	}
	public Point findTSS(){
		if(strand=='+'){
			TUnit fExon = exons.first();
			TSS = new Point(fExon.getCoords().getGenome(), fExon.getCoords().getChrom(), fExon.getCoords().getStart());
		}else{
			TUnit fExon = exons.last();
			TSS = new Point(fExon.getCoords().getGenome(), fExon.getCoords().getChrom(), fExon.getCoords().getEnd());
		}
		return(TSS);
	}
	
	
	public int getExonLength(){
		int eLen = 0;
		for(TUnit e : exons){
			eLen += e.getCoords().getWidth();
		}
		return(eLen);
	}
	
	public double getExonHitCount(){
		double eHits=0;
		for(TUnit e : exons){
			eHits += e.getHits();
		}
		return(eHits);
	}
	
	public boolean equals(Object o) { 
        if(!(o instanceof SplicedTUnit)) { return false; }
        SplicedTUnit g = (SplicedTUnit)o;
        if(!super.equals(g)) { return false; }
        if(exons != null || g.exons != null) { 
            if(exons != null && g.exons != null) { 
                if(exons.size() != g.exons.size()) { return false; }
                for(TUnit exon : exons) { if(!g.exons.contains(exon)) { return false; } }
            } else { 
                return false;
            }
        }
        return true;
    }
}
