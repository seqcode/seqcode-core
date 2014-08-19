package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;




/**
 * Similar to bedtools intersect with several added functions and access to readdb experiments
 * @author akshaykakumanu
 *
 */
public class Intersect {

	List<Point> a;
	List<Point> b;
	
	List<Map<Integer,Point>> aLocs;
	List<Map<Integer,Point>> bLocs;
	
	Integer[][] asort;
	Integer[][] bsort;
	
	Map<String, Integer> chrom2Id;
	Map<Integer, String> Id2Chrom;
	
	SeqLocator expta;
	SeqLocator exptb;
	
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
