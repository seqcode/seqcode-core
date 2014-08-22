package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.File;
import java.io.FileWriter;
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
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;




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
	
	private List<int[]> asort;
	private List<int[]> bsort;
	
	private Map<String, Integer> chrom2Id = new HashMap<String, Integer>();
	private Map<Integer, String> Id2Chrom = new HashMap<Integer,String>();

	private List<int[][]> aCounts;
	private List<int[][]> bCounts;
	
	/**
	 * replicate names are hard-coded at the moment, until I change this to the standard expt loading procedure
	 * replicates are for Zaret Mnase project
	 */
	private List<String> repsa = new ArrayList<String>();
	private List<String> repsb = new ArrayList<String>();
	
	Genome gen;
	
	int min_match_distance;
	
	public Intersect(List<Point> apeaks, List<Point> bpeaks, int min_match, Genome g) {
		a = apeaks;
		b= bpeaks;
		gen = g;
		
		int numChrom = 0;
		
		for(String chr : gen.getChromList()){
			if(!chr.contains("random")){
				chrom2Id.put(chr, numChrom);
				Id2Chrom.put(numChrom, chr);
				numChrom++;
			}
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
		
		asort = new ArrayList<int[]>();
		bsort = new ArrayList<int[]>();
		
		
		List<List<Integer>> tmp = new ArrayList<List<Integer>>();
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			tmp.add(new ArrayList<Integer>());
		}
		
		for( Point p : a ){
			tmp.get(chrom2Id.get(p.getChrom())).add(p.getLocation());
		}
		
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			int[] tmplocs = new int[tmp.get(c).size()];
			for(int l=0; l<tmp.get(c).size(); l++){
				tmplocs[l] = tmp.get(c).get(l);
			}
			int[] inds = StatUtil.findSort(tmplocs);
			tmplocs = StatUtil.permute(tmplocs,inds);
			asort.add(tmplocs);
		}
		
		tmp = null;
		
		tmp = new ArrayList<List<Integer>>();
		
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			tmp.add(new ArrayList<Integer>());
		}
		
		for( Point p : b ){
			tmp.get(chrom2Id.get(p.getChrom())).add(p.getLocation());
		}
		
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			int[] tmplocs = new int[tmp.get(c).size()];
			for(int l=0; l<tmp.get(c).size(); l++){
				tmplocs[l] = tmp.get(c).get(l);
			}
			int[] inds = StatUtil.findSort(tmplocs);
			tmplocs = StatUtil.permute(tmplocs,inds);
			bsort.add(tmplocs);
		}
		
		this.min_match_distance = min_match;
		
		
	}
	
	
	
	public Intersect(List<Point> apeaks, List<Point> bpeaks, Genome g, int minFragLen, int maxFragLen, int win, String expta, String exptb) throws SQLException, IOException {
		a = apeaks;
		b = bpeaks;
		gen = g;
		
		Integer numChrom = 0;
		
		for(String chr : gen.getChromList()){
			if(!chr.contains("random")){
				chrom2Id.put(chr, numChrom);
				Id2Chrom.put(numChrom, chr);
				numChrom++;
			}
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
		
		asort = new ArrayList<int[]>();
		bsort = new ArrayList<int[]>();
		
		
		List<List<Integer>> tmp = new ArrayList<List<Integer>>();
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			tmp.add(new ArrayList<Integer>());
		}
		
		for( Point p : a ){
			tmp.get(chrom2Id.get(p.getChrom())).add(p.getLocation());
		}
		
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			int[] tmplocs = new int[tmp.get(c).size()];
			for(int l=0; l<tmp.get(c).size(); l++){
				tmplocs[l] = tmp.get(c).get(l);
			}
			int[] inds = StatUtil.findSort(tmplocs);
			tmplocs = StatUtil.permute(tmplocs,inds);
			asort.add(tmplocs);
		}
		
		tmp = null;
		
		tmp = new ArrayList<List<Integer>>();
		
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			tmp.add(new ArrayList<Integer>());
		}
		
		for( Point p : b ){
			tmp.get(chrom2Id.get(p.getChrom())).add(p.getLocation());
		}
		
		for(int c=0; c<chrom2Id.keySet().size(); c++){
			int[] tmplocs = new int[tmp.get(c).size()];
			for(int l=0; l<tmp.get(c).size(); l++){
				tmplocs[l] = tmp.get(c).get(l);
			}
			int[] inds = StatUtil.findSort(tmplocs);
			tmplocs = StatUtil.permute(tmplocs,inds);
			bsort.add(tmplocs);
		}
		
		
		repsa.add("1a"); repsa.add("1b"); repsa.add("1c"); repsa.add("2a"); repsa.add("2b"); repsa.add("2c");
		repsb.add("1a"); repsb.add("1b"); repsb.add("1c"); repsb.add("2a"); repsb.add("2b"); repsb.add("2c");
		
		this.fillRepCounts(expta, exptb, minFragLen, maxFragLen, win);
		
	}
	
	private void fillRepCounts(String expta, String exptb, int minFragLen, int maxFragLen, int win) throws SQLException, IOException{
		String[] anames = expta.split(";");
		String[] bnames = exptb.split(";");
		System.out.println("Filling a counts");
		this.aCounts = this.getRepWindowCounts(anames[0], anames[1], repsa, asort, aLocs, minFragLen, maxFragLen, win);
		System.out.println(this.aCounts.get(chrom2Id.get("1")).length);
		System.out.println(this.aCounts.get(chrom2Id.get("3")).length);
		System.out.println("Filling b counts");
		this.bCounts = this.getRepWindowCounts(bnames[0], bnames[1], repsb, bsort, bLocs, minFragLen, maxFragLen, win);
	}
	
	
	
	private List<int[][]> getRepWindowCounts(String ename, String aname, List<String> reps, List<int[]> pointsSort, List<Map<Integer,Point>> pointsMap, int minFragLen, int maxFragLen, int win) throws SQLException, IOException{
		
		List<int[][]> ret = new ArrayList<int[][]>();
		
		for(int c=0; c< chrom2Id.size(); c++){
			
			ret.add(new int[pointsSort.get(c).length][reps.size()]);
			
			for(int r=0; r< reps.size(); r++){
				SeqLocator tmpLoc = new SeqLocator(ename, reps.get(r), aname);
				SeqExpander tmpExp = new SeqExpander(tmpLoc);
				Genome.ChromosomeInfo s = gen.getChrom(Id2Chrom.get(c));
				NamedRegion chrom = new NamedRegion(gen,Id2Chrom.get(c),1,s.getLength(),Id2Chrom.get(c));
				//System.out.println(this.gen.getChromID(Id2Chrom.get(c))+"\t"+Id2Chrom.get(c));
				Iterator<SeqHitPair> itr = tmpExp.getPairs(chrom);
				List<Integer> midpts = new ArrayList<Integer>();
				while(itr.hasNext()){
					SeqHitPair tmp = itr.next();
					if(tmp.getCode() == 1 && tmp.getMidpoint() != null){
						int  fragLen = Math.abs(tmp.getRight().getFivePrime() - tmp.getLeft().getFivePrime())+1;
						if(fragLen >= minFragLen && fragLen < maxFragLen){
							midpts.add(tmp.getMidpoint().getLocation());
						}
					}
					
				}
				
				int[] midsort = new int[midpts.size()];
				
				
				for(int i=0; i<midpts.size(); i++){
					midsort[i] = midpts.get(i);
				}
				
				int[] ind = StatUtil.findSort(midsort);
				midsort = StatUtil.permute(midsort, ind);
				
				
				for(int p=0; p< pointsSort.get(c).length; p++){
					ret.get(c)[p][r] = this.countRegionInChrom(pointsMap.get(c).get(pointsSort.get(c)[p]), midsort, win);
				}
				
				tmpExp.close();
			}
		}
		
		return ret;
	}
	
	private int countRegionInChrom(Point pt, int[] midsort, int win){
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
	 * @throws IOException 
	 */
	
	
	public void printWindCounts(String outbase) throws IOException{
		String dir = System.getProperty("user.dir");
		StringBuilder sb = new StringBuilder();
		FileWriter afw = new FileWriter(sb.append(dir).append("/").append(outbase).append("_a_Rep_counts.tab").toString());
		for(int c=0; c< aCounts.size(); c++){
			
			for(int p=0; p<aCounts.get(c).length; p++){
				StringBuilder o = new StringBuilder();
				o.append(aLocs.get(c).get(asort.get(c)[p]).getLocationString()).append("\t");
				for(int r=0; r<repsa.size(); r++){
					o.append(aCounts.get(c)[p][r]).append("\t");
				}
				afw.write(o.append("\n").toString());
			}
			
		}
		
		afw.close();
		StringBuilder sbb = new StringBuilder();
		FileWriter bfw = new FileWriter(sbb.append(dir).append("/").append(outbase).append("_b_Rep_counts.tab").toString());
		for(int c=0; c<bCounts.size(); c++){
			
			for(int p=0; p < bCounts.get(c).length; p++){
				StringBuilder ob = new StringBuilder();
				ob.append(bLocs.get(c).get(bsort.get(c)[p]).getLocationString()).append("\t");
				for(int r=0; r<repsb.size(); r++){
					ob.append(bCounts.get(c)[p][r]).append("\t");
				}
				bfw.write(ob.append("\n").toString());
			}
		}
		bfw.close();
	}
	
	
	public void intersect(){
		for(int c = 0; c<chrom2Id.keySet().size(); c++){
			for(int p = 0; p<asort.get(c).length; p++){
				int tmpKey = asort.get(c)[p];
				int match = Arrays.binarySearch(bsort.get(c), tmpKey);
				if( match < 0 ) { match = -match - 1; }
				int dist = 0;
				int bkey=0;
				if(match !=0){
					dist = ((bsort.get(c)[match]-tmpKey) < (tmpKey -bsort.get(c)[match-1]) )? bsort.get(c)[match]-tmpKey :tmpKey -bsort.get(c)[match-1];
					bkey = ((bsort.get(c)[match]-tmpKey) < (tmpKey -bsort.get(c)[match-1]) ) ? bsort.get(c)[match] : bsort.get(c)[match-1];
				}else{
					dist = 0;
					bkey = bsort.get(c)[match];
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
	
	
	public static void main(String[] args) throws NotFoundException, SQLException, IOException{
		ArgParser ap = new ArgParser(args);
		Genome g = Args.parseGenome(args).cdr();
		
		String aexpt = ap.getKeyValue("aExpt");
		String bexpt = ap.getKeyValue("bExpt");
		String apeaks;
		String bpeaks;
		int win;
		
		apeaks = ap.getKeyValue("aPeaks");
		win = Args.parseInteger(args, "win", -1);
		bpeaks = ap.getKeyValue("bPeaks");
			
		
		
		int minFragLen = Args.parseInteger(args, "minFrag", 140);
		int maxFragLen = Args.parseInteger(args, "maxFrag", 200);
		
		List<Point> apoints = Utilities.loadPeaksFromPeakFile(g, apeaks, win);
		List<Point> bpoints = Utilities.loadPeaksFromPeakFile(g, bpeaks, win);
		
		//Intersect analyzer = new Intersect(apoints, bpoints, g, minFragLen, maxFragLen, win, aexpt , bexpt);
		if(ap.hasKey("interesect")){
			
			int min_distance = Args.parseInteger(args, "minD", 40);
			Intersect analyzer = new Intersect(apoints, bpoints, min_distance, g);
			analyzer.intersect();
		}
		if(ap.hasKey("count")){
			Intersect analyzer = new Intersect(apoints, bpoints, g, minFragLen, maxFragLen, win, aexpt , bexpt);
			String outbase = ap.getKeyValue("out");
			analyzer.printWindCounts(outbase);
		}
		
	}
	
}
