package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;

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
	            chr = chr.replaceFirst("^chr", "");
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
			fivePrimeMethPerc = new float[numChroms][2][];
			
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
					this.fivePrimeMethPerc[chrom2ID.get(chr)][0]=null;
					this.fivePrimeMethPerc[chrom2ID.get(chr)][1]=null;
				}
			}
			
			//Sorting the arrays
			for(int i = 0; i < fivePrimePos.length; i++) {  // chr
				for(int j = 0; j < fivePrimePos[i].length; j++){
					if(fivePrimePos[i][j]!=null && fivePrimeCount[i][j]!=null && this.fivePrimeMethPerc[i][j] != null ){
						int[] inds = StatUtil.findSort(fivePrimePos[i][j]);
						fivePrimeCount[i][j] = StatUtil.permute(fivePrimeCount[i][j], inds);
						this.fivePrimeMethPerc[i][j] = StatUtil.permute(this.fivePrimeMethPerc[i][j], inds);
					}
				}
			}
			
			//free memenory by cleaning all the lists
			
			this.fivePrimePosList = null;
			this.fivePrimeCountList = null;
			this.fivePrimeMethPercList = null;
		
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
	
	private Pair<Double,Integer> getStrandedMethPerc(Region r, char strand){
		double totalMeth = 0;
		
		int num_sites = 0;
		String chr = r.getChrom();
		int chrID = chrom2ID.get(chr);
		int j = (strand == '+') ? 0:1;
		if(this.fivePrimePos[chrID][j] != null){
			int[] tempStarts = this.fivePrimePos[chrID][j];
			if(tempStarts.length !=0){
				int start_ind = Arrays.binarySearch(tempStarts, r.getStart());
				int end_ind   = Arrays.binarySearch(tempStarts, r.getEnd());
				if( start_ind < 0 ) { start_ind = -start_ind - 1; }
				if( end_ind < 0 )   { end_ind   = -end_ind - 1; }
				
				while (start_ind > 0 && tempStarts[start_ind - 1] >= r.getStart() ) {
	                start_ind--;
	            }
				while (end_ind < tempStarts.length && tempStarts[end_ind] <= r.getEnd()) {
	                end_ind++;
	            }
				
				for(int k = start_ind; k < end_ind; k++) {
					totalMeth += this.fivePrimeMethPerc[chrID][j][k];
					num_sites++;
				}
			}
		}
		
		
		return new Pair((num_sites == 0)? 0.0 : totalMeth/num_sites , num_sites);
	}
	
	public void initialize(){
		this.fivePrimePosList = new HashMap<String, ArrayList<Integer>[]>();
		this.fivePrimeCountList = new HashMap<String, ArrayList<Integer>[]>();
		this.fivePrimeMethPercList = new HashMap<String, ArrayList<Float>[]>();
	}
	
	
	//gettors
	
	public Pair<Double, Integer> getMethPerc(Region r){
		double ret=0;
		double pos_meth = this.getStrandedMethPerc(r, '+').car();
		double neg_meth = this.getStrandedMethPerc(r, '-').car();
		
		double tot_meth = pos_meth*this.getStrandedMethPerc(r, '+').cdr() + neg_meth*this.getStrandedMethPerc(r, '-').cdr();
		int tot_sites = this.getStrandedMethPerc(r, '+').cdr() + this.getStrandedMethPerc(r, '-').cdr();
		
		return new Pair((tot_sites == 0 ? 0 : tot_meth/tot_sites),tot_sites);
		
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
	
	

}
