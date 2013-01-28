package edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;

/*
 * GeneTUnit defines a gene that contains one to many transcripts.
 */
public class GeneTUnit extends SplicedTUnit{

	protected HashMap<String, SplicedTUnit> isoforms = new HashMap<String, SplicedTUnit>();
	protected TreeSet<Point> TSSs = new TreeSet<Point>();
	
	public GeneTUnit(String name, String ID, String type, Region r, char str) {
		super(name, ID, type, r, str);
	}
	
	//Accessors
	public int getNumIsoforms(){return isoforms.size();}
	public Collection<SplicedTUnit> getIsoforms(){return isoforms.values();}
	public Collection<Point> getTSSs(){return TSSs.descendingSet();}
	
	public void addIsoform(SplicedTUnit i){
		isoforms.put(i.ID, i);
		
		Iterator<TUnit> ite = i.getExonIterator();
		while(ite.hasNext()){
			this.addExon(ite.next());
		}
		
		Iterator<TUnit> itc = i.getCDSIterator();
		while(itc.hasNext()){
			this.addCDS(itc.next());
		}
		this.findIntrons();
		
		Point currTSS = i.getTSS();
		if(currTSS!=null && !TSSs.contains(currTSS)){
			TSSs.add(currTSS);
		}
	}
	public void addExon(String tID, TUnit c){
		if(hasIsoform(tID)){
			this.addExon(c);
			isoforms.get(tID).addExon(c);
		}
	}
	public void addCDS(String tID, TUnit c){
		if(hasIsoform(tID)){
			this.addCDS(c);
			isoforms.get(tID).addCDS(c);
		}
	}

	public boolean hasIsoform(String id){
		return isoforms.containsKey(id);
	}
	
	public SplicedTUnit getIsoform(String id){
		if(hasIsoform(id))
			return isoforms.get(id);
		else
			return null;
	}	
}
