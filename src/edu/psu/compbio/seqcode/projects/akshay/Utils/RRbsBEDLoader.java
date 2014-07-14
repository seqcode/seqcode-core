package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class RRbsBEDLoader {
	
	private HashMap<String, ArrayList<Integer>[]> fivePrimePosList = null;
	
	private HashMap<String, ArrayList<Integer>[]> fivePrimeCountList = null;
	
	private HashMap<String, ArrayList<Float>[]> fivePrimeMethPercList = null;
	
	private int[][][] fivePrimePos;
	private int[][][] fivePrimeCount;
	private float[][][] fivePrimeMethPerc;
	
	private HashMap<String, Integer> chrom2ID=new HashMap<String,Integer>();
	private HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
	
	private Genome gen;
	private File bedFile;
	
	
	private int totalReads;
	private int totalPos;
	
	public RRbsBEDLoader(Genome g, File f) {
		this.gen = g;
		this.bedFile = f;
	}
	
	public void sourceReads(){
		this.initialize();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(this.bedFile));
			String line;
			while((line = reader.readLine()) != null){
				line = line.trim();
				String[] words = line.split("\\s+");
				String chr="."; String strand = "";
	            int start=0, end=0;
	            chr = words[0];
	            start = Integer.parseInt(words[1]);
	            end = Integer.parseInt(words[2]);
	            strand = words[5];
	            int cov = Integer.parseInt(words[9]);
	            float percMeth = Float.parseFloat(words[10]);
	            if(!this.fivePrimePosList.containsKey(chr)){
	            	this.addChr(chr);
	            }
	            this.addHit(strand, start, end, chr, cov, percMeth);
	            this.totalReads = this.totalReads  +cov;
	            this.totalPos ++;
			}
			
			// Closing the file parser
			reader.close();
			
			int numChroms=0; 
			for(String chr : gen.getChromList()){
				chrom2ID.put(chr, numChroms);
				id2Chrom.put(numChroms, chr);
				numChroms++;
			}
			
			
			fivePrimePos  = new int[numChroms][2][];
			fivePrimeCount = new int[numChroms][2][];
			
			for(String chr : gen.getChromList()){
				if(this.fivePrimePosList.containsKey(chr)){
					for(int j = 0; j < fivePrimePosList.get(chr).length; j++)
						fivePrimePos[chrom2ID.get(chr)][j] = list2int(fivePrimePosList.get(chr)[j]);
				}else{
					fivePrimePos[chrom2ID.get(chr)][0]=null;
					fivePrimePos[chrom2ID.get(chr)][1]=null;
				}
			}
			for(String chr : gen.getChromList()){
				if(this.fivePrimeCountList.containsKey(chr)){
					for(int j = 0; j < this.fivePrimeCountList.get(chr).length; j++)
						this.fivePrimeCount[chrom2ID.get(chr)][j] = list2int(this.fivePrimeCountList.get(chr)[j]);
				}else{
					this.fivePrimeCount[chrom2ID.get(chr)][0]=null;
					this.fivePrimeCount[chrom2ID.get(chr)][1]=null;
				}
			}
			
			for(String chr : gen.getChromList()){
				if(this.fivePrimeMethPercList.containsKey(chr)){
					for(int j = 0; j < this.fivePrimeMethPercList.get(chr).length; j++)
						this.fivePrimeMethPerc[chrom2ID.get(chr)][j] = list2float(this.fivePrimeMethPercList.get(chr)[j]);
				}else{
					this.fivePrimeCount[chrom2ID.get(chr)][0]=null;
					this.fivePrimeCount[chrom2ID.get(chr)][1]=null;
				}
			}
			
			
			
			
			
			
			
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	private void addHit(String strand, int start, int end, String chr, int cov, float percMeth ){
		int strandID = strand.equals("+") ? 0 : 1;
		this.fivePrimePosList.get(chr)[strandID].add(strand.equals("+") ? start : end);
		this.fivePrimeCountList.get(chr)[strandID].add(cov);
		this.fivePrimeMethPercList.get(chr)[strandID].add(percMeth);
	}
	
	private void addChr(String chr){
		ArrayList<Integer>[] currIArrayList = new ArrayList[2];
		currIArrayList[0]=new ArrayList<Integer>();
		currIArrayList[1]=new ArrayList<Integer>();
		fivePrimePosList.put(chr, currIArrayList);
		ArrayList<Integer>[] currFArrayList = new ArrayList[2];
		currFArrayList[0]=new ArrayList<Integer>();
		currFArrayList[1]=new ArrayList<Integer>();
		fivePrimeCountList.put(chr, currFArrayList);
		ArrayList<Float>[] currFpArrayList = new ArrayList[2];
		currFpArrayList[0]=new ArrayList<Float>();
		currFpArrayList[1]=new ArrayList<Float>();
		this.fivePrimeMethPercList.put(chr, currFpArrayList);
	}
	
	public void initialize(){
		this.fivePrimePosList = new HashMap<String, ArrayList<Integer>[]>();
		this.fivePrimeCountList = new HashMap<String, ArrayList<Integer>[]>();
		this.fivePrimeMethPercList = new HashMap<String, ArrayList<Float>[]>();
	}
	
	
	//gettors
	
	public float getMethPerc(Region r){
		float ret=0;
		
		
		
		return ret;
	}
	
	protected float[] list2float(List<Float> list) {
		float[] out = new float[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	
	protected int[] list2int(List<Integer> list) {
		int[] out = new int[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	
	public HashMap<String, ArrayList<Integer>[]> getFivePos(){
		return this.fivePrimePosList;
	}
	
	public HashMap<String, ArrayList<Integer>[]> getFiveCounts(){
		return this.fivePrimeCountList;
	}
	public HashMap<String, ArrayList<Float>[]> getFiveMethPerc(){
		return this.fivePrimeMethPercList;
	}

}
