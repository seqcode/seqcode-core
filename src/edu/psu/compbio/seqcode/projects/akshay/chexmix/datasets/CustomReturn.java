package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

public class CustomReturn {
	public double pcc;
	public Vec maxVec1;
	public Vec maxVec2;
	
	public Vec maxvec;
	
	public BindingLocation bl;
	
	public int[] profile;
	public int no_in_seed;
	
	/**
	 * This constructor is used as the return object for scanBlWithBl method in the Bindinglocation class
	 * @param pcc
	 * @param maxVec1
	 * @param maxVec2
	 */
	public CustomReturn(double pcc, Vec maxVec1, Vec maxVec2) {
		this.pcc = pcc;
		this.maxVec1 = maxVec1;
		this.maxVec2 = maxVec2;
	}
	/**
	 * This constructor is used as the return object for ScanVecWithBl method in the BindingLocation class
	 * @param pcc
	 * @param maxvec
	 */
	public CustomReturn(double pcc, Vec maxvec){
		this.pcc = pcc;
		this.maxvec = maxvec;
	}
	/**
	 * This constructor is usesa as the return object in locationscanner class. The elements of the allblscan are these objects 
	 * @param pcc
	 * @param maxvec
	 * @param bl
	 */
	public CustomReturn(double pcc, Vec maxvec, BindingLocation bl){
		this.pcc = pcc;
		this.maxvec = maxvec;
		this.bl = bl;
	}
	
	/**
	 * used as the return object by the executescheme4 method of the buildseed class
	 * @param profile
	 * @param no_in_seed
	 */
	public CustomReturn(int[] profile, int no_in_seed){
		this.profile = profile;
		this.no_in_seed = no_in_seed;
	}
	

}
