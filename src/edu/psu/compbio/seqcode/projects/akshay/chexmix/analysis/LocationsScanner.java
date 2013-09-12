package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Config;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.CustomReturn;

public class LocationsScanner {
	
	public List<CustomReturn> scanOut = new ArrayList<CustomReturn>();
	public List<CustomReturn> allblscan = new ArrayList<CustomReturn>();
	private List<BindingLocation> failscanbls =  new ArrayList<BindingLocation>();
	private List<BindingLocation> passscanbls = new ArrayList<BindingLocation>();
	
	public LocationsScanner(List<BindingLocation> allbls, Config conf, int[] seedprofile) {
		for(BindingLocation bl : allbls){
			CustomReturn temp = bl.scanConcVecWithBl(seedprofile, conf.getIntSize());
			CustomReturn pushed = new CustomReturn(temp.pcc, temp.maxvec, bl);
			allblscan.add(pushed);
			if(temp.pcc > conf.getPccCutoff()){
				scanOut.add(pushed);
				passscanbls.add(bl);
			}
			else{
				failscanbls.add(bl);
			}
		}
	}
	
	//Accessorie
	
	public double[] getListofPCCvaluesThatPassCutoff(){
		double[] ret = new double[scanOut.size()];
		for(int i=0; i< scanOut.size(); i++){
			ret[i] = scanOut.get(i).pcc;
		}
		return ret;
	}
	
	public int[][] getTagsThatPassCuttoff(Config conf){
		int[][] ret = new int[scanOut.size()][4*conf.getIntSize()];
		for(int i=0; i<scanOut.size(); i++){
			List<Integer> addToRet = scanOut.get(i).bl.getConcatenatedTags(scanOut.get(i).maxvec.midpoint, conf.getIntSize(), scanOut.get(i).maxvec.orientation);
			for(int j=0; j<addToRet.size(); j++){
				ret[i][j] = addToRet.get(j);
			}
		}
		return ret;
	}
	
	public double[] getEntireListOfPccValues(){
		double[] ret =  new double[allblscan.size()];
		for(int i=0; i<allblscan.size(); i++){
			ret[i] = allblscan.get(i).pcc;
		}
		return ret;
	}
	
	public String[] getBlnamesThatPassCutoff(){
		String[] ret = new String[scanOut.size()];
		for(int i=0; i<scanOut.size(); i++){
			ret[i] = scanOut.get(i).bl.getName();
		}
		return ret;
	}
	
	public List<BindingLocation> getListOfBlsThatDoNotPassCuttOff(){
		return failscanbls;
	}
	
	public List<BindingLocation> getListOfBlsThatPassCuttoff(){
		return passscanbls;
	}
	
	public List<String> getNamesAndOrientationOfBlsThatPassCuttoff(){
		List<String> ret = new ArrayList<String>();
		for(CustomReturn cr: scanOut){
			String temp = cr.bl.getName()+":"+cr.maxvec.orientation;
			ret.add(temp);
		}
		return ret;
	}
	

}
