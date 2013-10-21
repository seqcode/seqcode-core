package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis.LoadTags;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.utils.*;

public class BindingLocation {

	/**
	 * seqpos and seqneg are Seq  (edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Seq) objects that contain the "+" and the "-" strand sequences 
	 * for the given bindinglocation 
	 */
	public Seq seqpos;
	public Seq seqneg;
	/**
	 * vecpos and vecneg are Vec objects (edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Vec) that contain the tags/reads in the "+" and "-" strands
	 */
	public Vec vecpos;
	public Vec vecneg;
	/**
	 * the chromosome (chr7, chr8 ..) in which the current bindinglocation falls
	 */
	private String chr;
	/**
	 * the midpoint (the actual genomic coordinate) of the current binding location
	 */
	private int midpoint;
	/**
	 * the range of the current binding location (the actual size of the location is 2*range) 
	 */
	private int range;
	/**
	 * the first element of this list is the start coordinate of the binding locaiton (inclusive)
	 * the second element of this list is the end coordinate of the binding location (not inclusive)
	 */
	private List<Integer> coords = new ArrayList<Integer>();
	
	private int Smoothsize;
	
	/**
	 * The only constructor for this class
	 * @param midpoint
	 * @param chr
	 * @param range
	 */
	public BindingLocation(int midpoint, String chr, Config conf) {
		this.chr = chr;
		this.midpoint = midpoint;
		this.range = conf.getBlsize();
		int start = midpoint - range;
		int end = midpoint + range;
		this.coords.add(start);
		this.coords.add(end);
		this.Smoothsize = conf.smoothing;
		
	}
	
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		BindingLocation bl = (BindingLocation) obj;
		return this.midpoint == bl.midpoint && this.chr == bl.chr && this.range == bl.range && this.Smoothsize == bl.Smoothsize;
	}
	
	@Override
	public int hashCode(){
		int result = 17;
		int code = (int) this.range;
		code+= (int) this.midpoint;
		code+= (int) this.Smoothsize;
		code += (int) (this.chr == null ? 0 :this.chr.hashCode());
		result = result*37 + code;
		return result;
	}
	
	public void fillSeqs(String posSeq){
		this.seqpos = new Seq(this.midpoint, this.range, this.chr, "+", posSeq);
		this.seqneg = new Seq(this.midpoint, this.range, this.chr, "-", SequenceUtils.reverseComplement(posSeq));
	}
	
	/**
	 * This method fills the vecpos and vecneg objects. Should be called only if you want to work with tags/tag distribution
	 * @param loader (LoadTags class object should be initiated first (it stores all the tags in a 3-d array format, see that class for more details))
	 */
	public void filltags(LoadTags loader){
		QueryTags tagsfetcher = new QueryTags(this.midpoint,this.range, this.chr, this.Smoothsize);
		//debug line
		//System.out.println(this.getName());
		//end
		this.vecpos = tagsfetcher.getTags(loader, "+");
		this.vecneg = tagsfetcher.getTags(loader, "-");
		if(this.Smoothsize >0){
			Smoothing smoo =  new Smoothing();
			Vec temppos = smoo.doSmoothing(this.vecpos, this.Smoothsize);
			if(temppos.tags.size() != vecpos.tags.size()){System.err.println("Somthing Wrong in the Smoonthing code; This error message is for the developer");}
			this.vecpos = temppos;
			smoo = new Smoothing();
			Vec tempneg = smoo.doSmoothing(this.vecneg, this.Smoothsize);
			if(tempneg.tags.size() != vecneg.tags.size()){System.err.println("Somthing Wrong in the Smoonthing code; This error message is for the developer");}
			this.vecneg = tempneg;
		}
	}
	
	/**
	 * Compares to the current binding location with the given binding location for all possible offsets and reversals and returns the offset and reversal that 
	 * gives the most highest similarity score for both binding locations.
	 * This method is only used while constructing the seed profile
	 * @param givenBL
	 * @param range
	 * @param smoothsize
	 * @return (Returns a custom object(edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.CustomReturn). maxvec1 in the return object is the offset
	 * and reversal that yields the maximum similarity in the current binding location. maxvec2 is for the given binding location)
	 */
	public CustomReturn scanBlWithBl(BindingLocation givenBL,int range){
		Vec maxVec1=null;
		Vec maxVec2=null;
		double pcc = -2.0;
		List<Integer> thisvec = this.getListMidpoints(range);
		List<Integer> givenvec = givenBL.getListMidpoints(range);

		
		for(int i=0; i<thisvec.size(); i++){
			for(int j=0; j<givenvec.size(); j++){
				List<Integer> first = this.getConcatenatedTags(thisvec.get(i), range, "+");
				List<Integer> second = givenBL.getConcatenatedTags(givenvec.get(j), range, "+");
				Pearson pccdriver = new Pearson(first,second);
				double temppcc = pccdriver.doComparision();
				if(temppcc>pcc ){
					pcc=temppcc;
					maxVec1= this.getSubVec(thisvec.get(i), range, "+");
					maxVec2= givenBL.getSubVec(givenvec.get(j), range, "+");
				}
				first = this.getConcatenatedTags(thisvec.get(i), range, "+");
				second = givenBL.getConcatenatedTags(givenvec.get(j),range, "-");
				pccdriver = new Pearson(first,second);
				temppcc = pccdriver.doComparision();
				if(temppcc>pcc ){
					//System.out.println(pcc);
					maxVec1= this.getSubVec(thisvec.get(i), range, "+");
					maxVec2= givenBL.getSubVec(givenvec.get(j), range, "-");
				}
				first = this.getConcatenatedTags(thisvec.get(i), range, "-");
				second = givenBL.getConcatenatedTags(givenvec.get(j), range, "+");
				pccdriver = new Pearson(first,second);
				temppcc = pccdriver.doComparision();
				if(temppcc>pcc ){
					//System.out.println(pcc);
					maxVec1= this.getSubVec(thisvec.get(i), range, "-");
					maxVec2= givenBL.getSubVec(givenvec.get(j), range, "+");
				}
				first = this.getConcatenatedTags(thisvec.get(i),range, "-");
				second = givenBL.getConcatenatedTags(givenvec.get(j), range, "-");
				pccdriver = new Pearson(first,second);
				temppcc = pccdriver.doComparision();
				if(temppcc>pcc ){
					pcc=temppcc;
					maxVec1= this.getSubVec(thisvec.get(i), range, "-");
					maxVec2= givenBL.getSubVec(givenvec.get(j), range, "-");
				}
			}
			
		}
		
		CustomReturn ret = new CustomReturn(pcc,maxVec1,maxVec2);
		return ret;
	}
	
	/**
	 * Given a vec (characterised by midpoint, range and orientation) in givenbl, this method finds a vector in the current BL that has the maximum similarity.
	 * This method is used by scheme 1 for building the seed profile
	 * @param givenbl
	 * @param midpoint
	 * @param orientation
	 * @param range
	 * @param smoothsize
	 * @return
	 */
	public CustomReturn scanVecWithBl(BindingLocation givenbl, int midpoint, String orientation, int range){
		List<Integer> first = givenbl.getConcatenatedTags(midpoint, range, orientation);
		List<Integer> thisvec = this.getListMidpoints(range);
		double maxpcc=-2.0;
		Vec maxvec = null;
		for(int i=0; i<thisvec.size(); i++){
			List<Integer> second = this.getConcatenatedTags(thisvec.get(i), range, "+");
			Pearson pcccalculater = new Pearson(first, second);
			double pcc = pcccalculater.doComparision();
			if(pcc>maxpcc){
				maxpcc = pcc;
				maxvec = this.getSubVec(thisvec.get(i), range, "+");
			}
			
			second = this.getConcatenatedTags(thisvec.get(i), range, "-");
			pcccalculater = new Pearson(first, second);
			pcc = pcccalculater.doComparision();
			if(pcc>maxpcc){
				maxpcc = pcc;
				maxvec = this.getSubVec(thisvec.get(i), range, "-");
			}
		}
		CustomReturn ret = new CustomReturn(maxpcc, maxvec);
		return ret;
	}
	
	/**
	 * 
	 * @param tags
	 * @param smoothsize
	 * @return
	 */
	public CustomReturn scanConcVecWithBl(int[] tags, int range){
		List<Integer> first = new ArrayList<Integer>();
		for(int j=0; j< tags.length; j++){
			first.add(tags[j]);
		}
		List<Integer> thisvec = this.getListMidpoints(range);
		double maxpcc=-2.0;
		Vec maxvec = null;
		for(int i=0; i<thisvec.size(); i++){
			List<Integer> second = this.getConcatenatedTags(thisvec.get(i), range, "+");
			Pearson pcccalculater = new Pearson(first, second);
			double pcc = pcccalculater.doComparision();
			if(pcc>maxpcc){
				maxpcc = pcc;
				maxvec = this.getSubVec(thisvec.get(i), range, "+");
			}
			
			second = this.getConcatenatedTags(thisvec.get(i), range, "-");
			pcccalculater = new Pearson(first, second);
			pcc = pcccalculater.doComparision();
			if(pcc>maxpcc){
				maxpcc = pcc;
				maxvec = this.getSubVec(thisvec.get(i), range, "-");
			}
			
		}
		CustomReturn ret = new CustomReturn(maxpcc, maxvec);
		return ret;
	}
	
	// Accessors
	
	/**
	 * Given a midpoint and range and orientaton, the methods returns the seq object if it is a subsequence of the current binding location
	 * @param midpoint
	 * @param orientation
	 * @param range
	 * @return
	 */
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
	
	/**
	 * Returns a concatenated integer list of tags for a given list of input paramenters
	 * @param midpoint
	 * @param chr
	 * @param range
	 * @param Orientation
	 * @param smoothsize
	 * @return
	 */
	public List<Integer> getConcatenatedTags(int midpoint, int range, String Orientation){
		List<Integer> ret  = new ArrayList<Integer>();
		if(Orientation  == "+"){
			Vec temp = this.getSubVec(midpoint, range, Orientation);
			Vec rev = this.getSubVec(midpoint, range, "-");
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
			Vec temp = this.getSubVec(midpoint, range, Orientation);
			Vec rev = this.getSubVec(midpoint, range, "+");
			List<Integer> tempvalues = new ArrayList<Integer>(temp.tags.values());
			List<Integer> revvalues = new ArrayList<Integer>(rev.tags.values());
			for(int i=tempvalues.size()-1; i>= 0; i--){
				ret.add(tempvalues.get(i));
			}
			for(int i=0;i<revvalues.size();i++){
				ret.add(revvalues.get(i));
			}
			
		}
		
		return ret;
	}
	
	/**
	 * Returns the vec object for a given input paraments.
	 * It is important to note that this method does not concatenate the tags from both strands. Use getConcatenatedTags to achieve that
	 * @param midpoint
	 * @param range
	 * @param orientation
	 * @param smoothsize
	 * @return
	 */
	public Vec getSubVec(int midpoint, int range, String orientation){
		Vec ret = null;
		if(orientation=="+"){
			ret = vecpos.getSub(midpoint, range);
		}
		if(orientation == "-"){
			ret =vecneg.getSub(midpoint, range);
		}
		
		return ret;
	}
	
	/**
	 * returns the name of the current binding location
	 * eg: chr8:777878, where 777878 is the midpoint and the chr7 is the chromosome
	 * @return
	 */
	public String getName(){
		return this.chr+":"+this.midpoint;
	}
	
	/**
	 * gives the list of vec objects that are a subset of the given binding location
	 * @param range
	 * @param orientation
	 * @param smoothsize
	 * @return
	 */
	public List<Vec> getListSubVec(int range, String orientation){
		List<Vec> ret = new ArrayList<Vec>();
		for(int i=this.coords.get(0)+range; i<this.coords.get(1)-range; i++){
			Vec temp = this.getSubVec(i, range, orientation);
			ret.add(temp);
		}
		return ret;
	}
	
	/**
	 * Similar to getListSubVec but just returns the list of midpoints.
	 * @param range
	 * @return
	 */
	public List<Integer> getListMidpoints(int range){
		List<Integer> ret = new ArrayList<Integer>();
		for(int i=this.coords.get(0)+range; i< this.coords.get(1)-range; i++){
			ret.add(i);
		}
		return ret;
	}
	
	public String getChr(){return this.chr;};
	public List<Integer> getCoords(){return this.coords;}
	
	
	
	
	
	
}
