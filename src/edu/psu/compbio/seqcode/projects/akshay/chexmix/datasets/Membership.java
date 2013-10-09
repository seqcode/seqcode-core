package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.util.List;

public class Membership {
	
	public BindingLocation bl;
	public CustomReturn cr;
	public int membership;
	
	public Membership(BindingLocation bl, CustomReturn cr, int membership) {
		this.bl = bl;
		this.cr = cr;
		this.membership = membership;
	}
	
	public double getPCC(){ return this.cr.pcc;}
	public int[] getConcVecAtMaxPCC(){
		List<Integer> temp = this.bl.getConcatenatedTags(this.cr.maxvec.midpoint, this.cr.maxvec.range, this.cr.maxvec.orientation);
		int[] ret = new int[temp.size()];
		for(int i=0; i< temp.size(); i++){
			ret[i] = temp.get(i);
		}
		return ret;
	}
	public String getPointName(){return this.cr.maxvec.chr+":"+this.cr.maxvec.midpoint+":"+this.cr.maxvec.orientation;}

}
