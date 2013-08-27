package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis.LoadTags;
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
	
	
	public BindingLocation(int midpoint, String chr, int range) {
		this.chr = chr;
		this.midpoint = midpoint;
		int start = midpoint - range;
		int end = midpoint + range;
		this.coords.add(start);
		this.coords.add(end);
		this.range = range;
	}
	
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		BindingLocation bl = (BindingLocation) obj;
		return this.midpoint == bl.midpoint && this.chr == bl.chr && this.range == bl.range;
	}
	
	@Override
	public int hashCode(){
		int result = 17;
		int code = (int) this.range;
		code+= (int) this.midpoint;
		code += (int) (this.chr == null ? 0 :this.chr.hashCode());
		result = result*37 + code;
		return result;
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
		if(genome == "mm10"){
			QueryMm10 seqloader = new QueryMm10(this.chr,this.midpoint,this.range);
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
	
	public void filltags(LoadTags loader) throws IOException{
		QueryTags tagsfetcher = new QueryTags(this.midpoint,this.range, this.chr);
		this.vecpos = tagsfetcher.getTags(loader, "+");
		this.vecneg = tagsfetcher.getTags(loader, "-");
	}
	
	public List<Integer> getConcatenatedTags(int midpoint, String chr, int range, String Orientation, int smoothsize){
		List<Integer> ret  = new ArrayList<Integer>();
		if(Orientation  == "+"){
			Vec temp = this.getSubVec(midpoint, range, Orientation, smoothsize);
			Vec rev = this.getSubVec(midpoint, range, "-", smoothsize);
			List<Integer> tempvalues = new ArrayList<Integer>(temp.tags.values());
			List<Integer> revvalues = new ArrayList<Integer>(rev.tags.values());
			for(int i =0; i<tempvalues.size(); i++){
				ret.add(tempvalues.get(i));
			}
			for(int i= revvalues.size()-1;i>=0 ;i--){
				ret.add(revvalues.get(i));
				
			}
		}
		
		if(Orientation == "-"){
			Vec temp = this.getSubVec(midpoint, range, Orientation, smoothsize);
			Vec rev = this.getSubVec(midpoint, range, "+", smoothsize);
			List<Integer> tempvalues = new ArrayList<Integer>(temp.tags.values());
			List<Integer> revvalues = new ArrayList<Integer>(rev.tags.values());
			for(int i=revvalues.size()-1;i>=0;i--){
				ret.add(revvalues.get(i));
			}
			for(int i=0; i< tempvalues.size(); i++){
				ret.add(tempvalues.get(i));
			}
		}
		
		return ret;
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
	
	public String getName(){
		return this.chr+":"+this.midpoint;
	}
	
	public List<Vec> getListSubVec(int range, String orientation, int smoothsize){
		List<Vec> ret = new ArrayList<Vec>();
		for(int i=this.coords.get(0)+range; i<this.coords.get(1)-range; i++){
			Vec temp = this.getSubVec(i, range, orientation, smoothsize);
			ret.add(temp);
		}
		return ret;
	}
	
	
	public List<Integer> getListMidpoints(int range){
		List<Integer> ret = new ArrayList<Integer>();
		for(int i=this.coords.get(0)+range; i< this.coords.get(1)-range; i++){
			ret.add(i);
		}
		return ret;
	}
	
	
	/* maxVec1 is from thisvector and mxvec2 is from given vector
	 * */
	public CustomReturn scanTwoBLs(BindingLocation givenBL,int range, int smoothsize){
		Vec maxVec1=null;
		Vec maxVec2=null;
		double pcc = -2.0;
		List<Integer> thisvec = this.getListMidpoints(range);
		List<Integer> givenvec = givenBL.getListMidpoints(range); 
		
		for(int i=0; i<thisvec.size(); i++){
			for(int j=0; j<givenvec.size(); j++){
				List<Integer> first = this.getConcatenatedTags(thisvec.get(i), this.chr, range, "+", smoothsize);
				List<Integer> second = givenBL.getConcatenatedTags(givenvec.get(j), givenBL.chr, range, "+", smoothsize);
				//Pearson pccdriver = new Pearson(first,second);
				//double temppcc = pccdriver.doComparision();
				//if(temppcc>pcc ){
				//	pcc=temppcc;
					// debuggining lines
				//	System.out.println(pcc);
				//	maxVec1= this.getSubVec(thisvec.get(i), range, "+", smoothsize);
				//	maxVec2= givenBL.getSubVec(givenvec.get(j), range, "+", smoothsize);
				//}
				first = this.getConcatenatedTags(thisvec.get(i), this.chr, range, "+", smoothsize);
				second = givenBL.getConcatenatedTags(givenvec.get(j), givenBL.chr, range, "-", smoothsize);
				//pccdriver = new Pearson(first,second);
				//temppcc = pccdriver.doComparision();
				//if(temppcc>pcc ){
				//	pcc=temppcc;
					// debuggining lines
				//	System.out.println(pcc);
				//	maxVec1= this.getSubVec(thisvec.get(i), range, "+", smoothsize);
				//	maxVec2= givenBL.getSubVec(givenvec.get(j), range, "-", smoothsize);
				//}
				first = this.getConcatenatedTags(thisvec.get(i), this.chr, range, "-", smoothsize);
				second = givenBL.getConcatenatedTags(givenvec.get(j), givenBL.chr, range, "+", smoothsize);
				//pccdriver = new Pearson(first,second);
				//temppcc = pccdriver.doComparision();
				//if(temppcc>pcc ){
				//	pcc=temppcc;
					// debuggining lines
				//	System.out.println(pcc);
				//	maxVec1= this.getSubVec(thisvec.get(i), range, "-", smoothsize);
				//	maxVec2= givenBL.getSubVec(givenvec.get(j), range, "+", smoothsize);
				//}
				first = this.getConcatenatedTags(thisvec.get(i), this.chr, range, "-", smoothsize);
				second = givenBL.getConcatenatedTags(givenvec.get(j), givenBL.chr, range, "-", smoothsize);
				//pccdriver = new Pearson(first,second);
				//temppcc = pccdriver.doComparision();
				//if(temppcc>pcc ){
				//	pcc=temppcc;
					// debuggining lines
				//	System.out.println(pcc);
				//	maxVec1= this.getSubVec(thisvec.get(i), range, "-", smoothsize);
				//	maxVec2= givenBL.getSubVec(givenvec.get(j), range, "-", smoothsize);
				//}
			}
		}
		
		CustomReturn ret = new CustomReturn(pcc,maxVec1,maxVec2);
		
 		
		return ret;
	}
	
	
	
	
	
}
