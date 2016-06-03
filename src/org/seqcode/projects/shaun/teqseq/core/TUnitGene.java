package org.seqcode.projects.shaun.teqseq.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.seqcode.gse.projects.readdb.PairedHit;


public class TUnitGene extends TUnit{
	protected List<TUnit> componentUnits = new ArrayList<TUnit>();
	
	public TUnitGene(TUnit u, int i){
		super(u.getCoords(), i, u.getStrand());
		initializeCounters(u.getRepCount());
		addComponent(u);
	}
	
	public List<TUnit> getComponents(){return componentUnits;}
	public void sortComponents(){ Collections.sort(componentUnits);}
	
	public void addComponent(TUnit c){
		componentUnits.add(c);
		//Add coverage & starts
		addCoverageAndStarts(c);
		//Expand coords
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
		//Add annotations
		if(c.getAnnotations()!=null){
			for(String a : c.getAnnotations()){
				if(!annotations.contains(a))
					annotations.add(a);
			}
		}
	}

	public List<AlignHit> getUnitBoundaryHits(){
		List<AlignHit> hits = new ArrayList<AlignHit>();
		for(TUnit u :componentUnits)
			if(u.getBoundaryHits()!=null && u.getBoundaryHits().size()>0)
				hits.addAll(u.getBoundaryHits());
		return hits;
	}
	
	public List<PairedHit> getUnitBoundaryPairedHits(){
		List<PairedHit> hits = new ArrayList<PairedHit>();
		for(TUnit u :componentUnits)
			if(u.getBoundaryPairedHits()!=null && u.getBoundaryPairedHits().size()>0)
				hits.addAll(u.getBoundaryPairedHits());
		return hits;
	}
	
	protected void addCoverageAndStarts(TUnit c){
		for(int r=0; r<c.getRepCount(); r++){
			addCoverage(r, c.getHitStarts(r));
			addHitStart(r, c.getHitStarts(r));
		}
	}
	
	public boolean overlaps(AlignHit hit){
		boolean oLap = false;
		for(TUnit tu : componentUnits)
			oLap = oLap | tu.overlaps(hit);
		return(oLap);
	}
	
	public TUnitGene joinGene(TUnitGene g){
		for(TUnit u : g.getComponents())
			this.addComponent(u);
		return this;
	}
}
