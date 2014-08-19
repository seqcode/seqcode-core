package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHitPair;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;




/**
 * Similar to bedtools intersect with several added functions and access to readdb experiments
 * @author akshaykakumanu
 *
 */
public class Intersect {

	private List<Point> a;
	private List<Point> b;
	
	private List<Map<Integer,Point>> aLocs;
	private List<Map<Integer,Point>> bLocs;
	
	private Integer[][] asort;
	private Integer[][] bsort;
	
	private Map<String, Integer> chrom2Id;
	private Map<Integer, String> Id2Chrom;
	
	//private SeqLocator expta;
	//private SeqLocator exptb;
	
	private Integer[][][] aCounts;
	private Integer[][][] bCounts;
	
	/**
	 * replicate names are hard-coded at the moment, until I change this to the standard expt loading procedure
	 * replicates are for Zaret Mnase project
	 */
	private Collection<String> repsa = new ArrayList<String>();
	private Collection<String> repsb = new ArrayList<String>();
	
	Genome gen;
	
	int min_match_distance;
	
	
	public Intersect(List<Point> apeaks, List<Point> bpeaks, Genome g) {
		a = apeaks;
		b = bpeaks;
		gen = g;
		
		Integer numChrom = 0;
		
		for(String chr : gen.getChromList()){
			chrom2Id.put(chr, numChrom);
			Id2Chrom.put(numChrom, chr);
			numChrom++;
		}
		
		aLocs = new ArrayList<Map<Integer,Point>>();
		bLocs = new ArrayList<Map<Integer, Point>>();
		
		for(String s : chrom2Id.keySet()){
			aLocs.add(new HashMap<Integer, Point>());
			bLocs.add(new HashMap<Integer, Point>());
		}
		
		for( Point p :  a ){
			aLocs.get(chrom2Id.get(p.getChrom())).put(p.getLocation(), p);
		}
		
		for( Point p : b ){
			bLocs.get(chrom2Id.get(p.getChrom())).put(p.getLocation(), p);
		}
		
		asort = new Integer[chrom2Id.keySet().size()][a.size()];
		bsort = new Integer[chrom2Id.keySet().size()][b.size()];
		
		int count = 0;
		for( Point p : a ){
			asort[chrom2Id.get(p.getChrom())][count] = p.getLocation();
			count++;
		}
		
		count = 0;
		for( Point p : b ){
			bsort[chrom2Id.get(p.getChrom())][count] = p.getLocation();
		}
		
		for(int i=0; i<chrom2Id.keySet().size(); i++){
			int[] indsa = StatUtil.findSort(asort[i]);
			asort[i] = StatUtil.permute(asort[i],indsa);
			
			int[] indsb =  StatUtil.findSort(bsort[i]);
			bsort[i] = StatUtil.permute(bsort[i],indsb);
		}
		
		repsa.add("1a"); repsa.add("1b"); repsa.add("1c"); repsa.add("2a"); repsa.add("2b"); repsa.add("2c");
		repsb.add("1a"); repsb.add("1b"); repsb.add("1c"); repsb.add("2a"); repsb.add("2b"); repsb.add("2c");
		
	}
	
	
	public void fillCounts(String Aename, String Aaname, String Bename, String Baname,  int minFragLen, int maxFragLen, int win) throws SQLException, IOException{
		
		this.aCounts = new Integer[chrom2Id.keySet().size()][a.size()][repsa.size()];
		this.bCounts = new Integer[chrom2Id.keySet().size()][b.size()][repsb.size()];
		
		for(int c=0; c< aLocs.size(); c++){
			
			for(int r=0; r< repsa.size(); r++){
				SeqLocator atmpLoc = new SeqLocator(Aename, repsa, Baname);
				SeqExpander atmpExp = new SeqExpander(atmpLoc);
				Genome.ChromosomeInfo s = gen.getChrom(Id2Chrom.get(c));
				NamedRegion chrom = new NamedRegion(gen,Id2Chrom.get(c),1,s.getLength(),Id2Chrom.get(c));
				Iterator<SeqHitPair> aitr = atmpExp.getPairs(chrom);
				List<Integer> amidpts = new ArrayList<Integer>();
				while(aitr.hasNext()){
					SeqHitPair tmp = aitr.next();
					if(tmp.getCode() == 1 && tmp.getMidpoint() != null){
						int  fragLen = Math.abs(tmp.getRight().getFivePrime() - tmp.getLeft().getFivePrime())+1;
						if(fragLen >= minFragLen && fragLen < maxFragLen){
							amidpts.add(tmp.getMidpoint().getLocation());
						}
					}
					
				}
				
				int[] amidsort = new int[amidpts.size()];
				
				
				for(int i=0; i<amidpts.size(); i++){
					amidsort[i] = amidpts.get(i);
				}
				for(int p=0; p< aLocs.get(c).keySet().size(); p++){
					aCounts[c][p][r] = this.countHits(aLocs.get(c).get(p), amidsort, win);
				}
				
				atmpExp.close();
				
				
				
				
			}
		}
	}
	
	public int countHits(Point pt, int[] midsort, int win){
		Region r = pt.expand(win);
		int ret_count;
		int start_match = Arrays.binarySearch(midsort, r.getStart());
		int end_match = Arrays.binarySearch(midsort, r.getEnd());
		if(start_match <0 ) {start_match = -start_match -1;}
		if(end_match <0) {end_match = -end_match-1;}
		
		while(start_match >0 && midsort[start_match-1] >= r.getStart()){
			start_match--;
		}
		while(end_match < midsort.length && midsort[end_match] <= r.getEnd()){
			end_match++;
		}
		 ret_count = end_match-start_match;
		 
		return ret_count;		                                                                  
	}
	
	
	
	/**
	 * Prints the matching overlapping win
	 */
	public void intersect(){
		for(int c = 0; c<chrom2Id.keySet().size(); c++){
			for(int p = 0; p<a.size(); p++){
				int tmpKey = asort[c][p];
				int match = Arrays.binarySearch(bsort[c], tmpKey);
				if( match < 0 ) { match = -match - 1; }
				int dist = 0;
				int bkey=0;
				if(match !=0){
					dist = ((bsort[c][match]-tmpKey) < (tmpKey -bsort[c][match-1]) )? bsort[c][match]-tmpKey :tmpKey -bsort[c][match-1];
					bkey = ((bsort[c][match]-tmpKey) < (tmpKey -bsort[c][match-1]) ) ? bsort[c][match] : bsort[c][match-1];
				}else{
					dist = 0;
					bkey = bsort[c][match];
				}
				
				if(dist<min_match_distance){
					System.out.println(aLocs.get(c).get(tmpKey).getLocation()+"\t"+
							bLocs.get(c).get(bkey).getLocation()+"\t"+Integer.toString(dist));
				}
				else{
					System.out.println(aLocs.get(c).get(tmpKey).getLocation()+"\t"+"Inf"+"\t"+"Inf");
				}
			}
		}
	}
	
}
