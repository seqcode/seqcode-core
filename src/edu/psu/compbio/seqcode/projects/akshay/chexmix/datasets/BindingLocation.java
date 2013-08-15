package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.utils.*;

public class BindingLocation {

	public Seq seqpos;
	public Seq seqneg;
	public Vec vecpos;
	public Vec vecneg;
	private String chr;
	private int midpoint;
	private int range;
	private List<Integer> coords = new ArrayList<Integer>();
	private Map<Integer,Integer> tags = new TreeMap<Integer,Integer>();
	
	public BindingLocation(int midpoint, String chr, int range) {
		this.chr = chr;
		this.midpoint = midpoint;
		int start = midpoint - range;
		int end = midpoint + range;
		this.coords.add(start);
		this.coords.add(end);
		this.range = range;
	}
	
	public void fillSeqs(String genome) throws IOException{
		if(genome == "hg19"){
			QueryHg19 seqloader = new QueryHg19(this.chr,this.midpoint, this.range);
			seqloader.fillGenomePath();
			this.seqpos = seqloader.getSeq("+");
			this.seqneg = seqloader.getSeq("-");
		}
		if(genome == "mm9"){
			QueryMm9 seqloader = new QueryMm9(this.chr,this.midpoint,this.range);
			seqloader.fillGenomePath();
			this.seqpos = seqloader.getSeq("+");
			this.seqneg = seqloader.getSeq("-");
		}
	}
	
	public Seq getSubSeq(int midpoint, String orientation, int range){
		Seq ret = null;
		if(orientation == "+"){
			ret = seqpos.getSub(midpoint, range);
		}
		else{
			ret = seqneg.getSub(midpoint, range);
		}
		return ret;
	}
	
	public void filltags(String tagsfiletype, String tagsfilepath) throws IOException{
		if(tagsfiletype == "bed"){
			QueryTagsBed tagsloader = new QueryTagsBed(this.midpoint, this.range, this.chr);
			tagsloader.prepareQureybed("+");
			tagsloader.fillTagsBedPath(tagsfilepath);
			this.vecpos = tagsloader.getTags("+");
			tagsloader.prepareQureybed("-");
			this.vecneg = tagsloader.getTags("-");
		}
		if(tagsfiletype == "idx"){
			QueryTagsIdx tagsloader = new QueryTagsIdx(this.midpoint,this.range,this.chr);
			tagsloader.prepareQureybed("+");
			tagsloader.fillTagsBedPath(tagsfilepath);
			this.vecpos = tagsloader.getTags("+");
			tagsloader.prepareQureybed("-");
			this.vecneg = tagsloader.getTags("-");
		}
	}
	
	public Vec getSubVec(int midpoint, int range, String orientation, int smoothsize){
		Vec ret = null;
		if(smoothsize > 0){
			if(orientation=="+"){
				Smoothing smoo = new Smoothing();
				Vec tempret = smoo.doSmoothing(this.vecpos, smoothsize);
				ret = tempret.getSub(midpoint, range);
			}
			if(orientation =="-"){
				Smoothing smoo = new Smoothing();
				Vec tempret = smoo.doSmoothing(this.vecneg, smoothsize);
				ret = tempret.getSub(midpoint, range);
			}
		}
		else{
			if(orientation=="+"){
				ret = vecpos.getSub(midpoint, range);
			}
			if(orientation == "-"){
				ret =vecneg.getSub(midpoint, range);
			}
		}
		return ret;
	}
	
	public List<Vec> getListSubVec(int range, String orientation, int smoothsize){
		List<Vec> ret = new ArrayList<Vec>();
		for(int i=this.coords.get(0)+range; i<this.coords.get(1)-range; i++){
			Vec temp = this.getSubVec(i, range, orientation, smoothsize);
			ret.add(temp);
		}
		return ret;
	}
	
	
	
}
