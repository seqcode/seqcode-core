package edu.psu.compbio.seqcode.projects.akshay.regulatorydomains;

public class GeneDomain {
	private String geneName;
	private int clusterMembership =0;
	private double foldChange;
	
	public GeneDomain(String gn, double fchange) {
		geneName =gn;
		foldChange = fchange;
	}
	
	//Settors
	public void setClusterInd(int index){clusterMembership = index;}
	
	
	// Gettors
	public int getClusterIndex(){return clusterMembership;}
	public String getGeneName(){return geneName;}
	public double getFoldChange(){return foldChange;}
	

}
