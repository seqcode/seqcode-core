package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.SeqBiasModel;

/**
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public abstract class ReadLoader {

	protected Genome gen;
	protected String conditionName="";
	protected String sourceName = ""; //File name or other descriptor
	protected String sourcePath="";   //File path or URL, etc
	protected Map<String, String> chromNameTranslator = new HashMap<String,String>();
	
	public ReadLoader(Genome g, Map<String, String> nameTrans){
		gen = g;
		if(gen==null)
			gen = this.estimateGenome();
		chromNameTranslator = nameTrans;
	}
	//Accessors
	public Genome getGenome(){return gen;}
	public String getSourceName(){return sourceName;}
	public String getConditionName(){return conditionName;}
	public String getSourcePath(){return sourcePath;}
	public void setNameTranslator(Map<String, String> nameTrans){chromNameTranslator = nameTrans;}
	
	//Get all hits overlapping a region 
	public abstract List<AlignHit> getOverlappingHits(Region r, SeqBiasModel bias);
	public List<AlignHit> getOverlappingHits(Region r){return(getOverlappingHits(r, null));}
	
	//Get all reads overlapping a region (from the stream)
	public abstract List<AlignHit> getOverlappingHitsFromStream(Region r, SeqBiasModel bias);
	public List<AlignHit> getOverlappingHitsFromStream(Region r){return(getOverlappingHitsFromStream(r, null));}
	
	//Estimate the genome from the data itself
	public abstract Genome estimateGenome();
	
	//Initialize an iterator that traverses all reads (actual iterator is specific to each ReadLoader)
	public abstract void initializeReadIterator();
	
	//Close the read iterator
	public abstract void closeReader();
	
	//Does the read iterator have another hit left?
	public abstract boolean iteratorHasNextHit();
	
	//Increment the iterator (return false if iterator has run out of reads)
	public abstract boolean nextRead();
	
	//Get the next AlignHit from the read iterator
	public abstract AlignHit getCurrHit();
	
	//Get the next read sequence from the read iterator
	public abstract byte[] getCurrReadSeq();
}
