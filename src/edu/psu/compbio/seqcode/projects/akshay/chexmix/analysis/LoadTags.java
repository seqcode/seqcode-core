package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.io.File;
import java.util.*;
import java.util.HashMap;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.projects.multigps.hitloaders.*;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;

public class LoadTags {
	
	private String tagsfilepath;
	
	private Genome gen; // this genome object is estimated for the data internally
	
	private int numChroms; 
	
	/** 
	 * Five prime ends of the read hits. <br>
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the coordinates of the hits
	 */
	public int[][][] fivePrimePos = null;
	
	private double totalHits; //totalHits is the sum of alignment weights
	
	/**
	 * Sum of read hit weights that corresponds to the 5' position
	 * First dimension represents the corresponding chromosome ID. <br>
	 * Second dimension represents the strand. 0 for '+', 1 for '-' <br>
	 * Third dimension contains the number of hits at corresponding start position 
	 */
	public float[][][] fivePrimeCounts = null;
	
	public HashMap<String, Integer> chrom2ID=new HashMap<String,Integer>();
	public HashMap<Integer,String> id2Chrom=new HashMap<Integer,String>();
	
	/**
	 * Constructor for this class,
	 * Takses the path of the tags file
	 */
	public LoadTags(String path) {
		this.tagsfilepath  =path;
	}
	
	
	/**
	 * Fills the tagsfile path field
	 * This method is called by the loadHits public method
	 * @param path get this value from the config
	 */
	
	private void fillTagsfilepathname(String path){
		this.tagsfilepath = path;
	}
	
	/**
	 * Gives the Hitloader for the given input taglsfile
	 * This method is called by the loadHits public method
	 * @param format Tags file format: get this from the config class
	 * @param useNonUnique Usually false. get this from config file
	 * @return The appropriate HitLoader
	 */
	private HitLoader getFileHitloader(String format, boolean useNonUnique){
		HitLoader currReader=null;
		File file = new File(this.tagsfilepath);
		if(!file.isFile()){System.err.println("File not found: "+file.getName());System.exit(1);}
		if(format.equals("SAM") || format.equals("BAM")){
			currReader = new SAMFileHitLoader(file,useNonUnique);
		}else if(format.equals("TOPSAM")){
			currReader = new TophatFileHitLoader(file,useNonUnique);
		}else if(format.equals("ELAND")){
			currReader = new ElandFileHitLoader(file,useNonUnique);
		}else if(format.equals("NOVO")){
			currReader = new NovoFileHitLoader(file,useNonUnique);
		}else if(format.equals("BOWTIE")){
			currReader = new BowtieFileHitLoader(file,useNonUnique);
		}else if(format.equals("BED")){
			currReader = new BEDFileHitLoader(file,useNonUnique);
		}else if(format.equals("IDX")){
			currReader = new IDXFileHitLoader(file,useNonUnique);
		}else{
		    System.err.println("Unknown file format: "+format);
		    System.exit(1);
		}
		return currReader;
	}
	
	public void loadHits(String path, String format, boolean useNonUnique ){
		this.fillTagsfilepathname(path);
		HitLoader hloader = this.getFileHitloader(format, useNonUnique);
		
		hloader.sourceReads();
		//Estimate the genome
		this.gen = this.estimateGenome(hloader.getFivePrimePositions());
		
		//Initialize chromosome name to id maps 
		// The chr string is like the following: "7", "8" for chromosomes 7 and 8
		this.numChroms=0;
		for(String chr : this.gen.getChromList()){
			this.chrom2ID.put(chr, numChroms);
			id2Chrom.put(numChroms, chr);
			numChroms++;
		}
		
		//Initialize the data structures
		this.fivePrimePos  = new int[this.numChroms][2][];
		this.fivePrimeCounts = new float[this.numChroms][2][];
		
		//Copy over the data
		for(String chr: this.gen.getChromList()){
			if(hloader.getFivePrimePositions().containsKey(chr)){
				for(int j =0 ; j< hloader.getFivePrimePositions().get(chr).length; j++){
					this.fivePrimePos[this.chrom2ID.get(chr)][j] = this.list2int(hloader.getFivePrimePositions().get(chr)[j]);
				}
			}
			else{
				this.fivePrimePos[this.chrom2ID.get(chr)][0]=null;
				this.fivePrimePos[this.chrom2ID.get(chr)][1]=null;
			}
		}
		
		for(String chr: this.gen.getChromList()){
			if(hloader.getFivePrimeCounts().containsKey(chr)){
				for(int j =0 ; j< hloader.getFivePrimeCounts().get(chr).length; j++){
					this.fivePrimeCounts[this.chrom2ID.get(chr)][j] = this.list2float(hloader.getFivePrimeCounts().get(chr)[j]);
				}
			}
			else{
				this.fivePrimeCounts[this.chrom2ID.get(chr)][0]=null;
				this.fivePrimeCounts[this.chrom2ID.get(chr)][1]=null;
			}
		}
		
		//Sort the arrays 
		for(int i = 0; i < this.fivePrimePos.length; i++) {  // chr
			for(int j = 0; j < this.fivePrimePos[i].length; j++) { // strand
				if(fivePrimePos[i][j]!=null && fivePrimeCounts[i][j]!=null){
					int[] inds = StatUtil.findSort(fivePrimePos[i][j]);
					fivePrimeCounts[i][j] = StatUtil.permute(fivePrimeCounts[i][j], inds);
				}
			}
		}
		
		//Collapse duplicate positions
		//Testing if the finished set of positions contains duplicates
		for(int i = 0; i < fivePrimePos.length; i++){
			for(int j = 0; j < fivePrimePos[i].length; j++){
				if(fivePrimePos[i][j]!=null && fivePrimePos[i][j].length>0){
					int uniquePos=1;
					for(int k = 0; k < fivePrimePos[i][j].length-1; k++)
						if(fivePrimePos[i][j][k+1]!=fivePrimePos[i][j][k]){uniquePos++;}
					int[] tmpPos = new int[uniquePos];
					float[] tmpCnt = new float[uniquePos];
					for(int x=0; x<uniquePos; x++){tmpCnt[x]=0;}
					int x=0;
					tmpPos[x] = fivePrimePos[i][j][0];
					tmpCnt[x] += fivePrimeCounts[i][j][0];
					for(int k = 1; k < fivePrimePos[i][j].length; k++){
						if(fivePrimePos[i][j][k]!=fivePrimePos[i][j][k-1]){x++;}
						tmpPos[x] = fivePrimePos[i][j][k];
						tmpCnt[x] += fivePrimeCounts[i][j][k];
					}
					fivePrimeCounts[i][j] = tmpCnt;
					fivePrimePos[i][j] = tmpPos;
				}
			}
		}
		
		this.updateTotalHits();
		hloader.resetLoader();
		
	}
	
	
	
	/**
	 * Estimate a genome from the observed read positions that are collected into the list
	 * @param posList HashMap indexed by chr containing ArrayLists of hit positions
	 * @return Genome
	 */
	private Genome estimateGenome(HashMap<String, ArrayList<Integer>[]> posList){
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		for(String c : posList.keySet()){
			int max = 0;
			for(int j=0; j<posList.get(c).length; j++){
				for(Integer i : posList.get(c)[j]){
					if(i>max)
						max=i;
				}
			}
			chrLenMap.put(c, max);
		}
		Genome g =new Genome("Genome", chrLenMap);
		return(g);
	}
	
	/**
	 * Recount hit weights
	 */
	private void updateTotalHits(){
		this.totalHits = 0.0;
		for(int i = 0; i < this.fivePrimeCounts.length; i++)
			for(int j = 0; j < this.fivePrimeCounts[i].length; j++)
				if(this.fivePrimeCounts[i][j]!=null)
					for(int k = 0; k < this.fivePrimeCounts[i][j].length; k++)
						this.totalHits += this.fivePrimeCounts[i][j][k];
	}
	
	
	/**
	 * Simple convertor of List to int[]
	 * @param list
	 * @return
	 */
	private int[] list2int(List<Integer> list) {
		int[] out = new int[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	
	/**
	 * Simple covertor of List to float[]
	 * @param list
	 * @return
	 */
	protected float[] list2float(List<Float> list) {
		float[] out = new float[list.size()];
		for(int i = 0; i < out.length; i++)
			out[i] = list.get(i);
		return out;
	}
	
	
	
	

}
