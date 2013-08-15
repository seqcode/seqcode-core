package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

public class QueryTagsBed extends QueryTags{

	public QueryTagsBed(int midpoint, int range, String chr) {
		super(midpoint, range, chr);
	}
	
	@Override
	public void fillTagsBedPath(String tagspath) {
		this.tagsbedfilepath = tagspath;
		
	}

}
