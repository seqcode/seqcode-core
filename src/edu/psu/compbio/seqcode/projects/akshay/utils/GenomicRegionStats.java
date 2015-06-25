package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class GenomicRegionStats {
	public List<List<Point>> points_lists = new ArrayList<List<Point>>();
	public List<String> factor_string = new ArrayList<String>();
	public int min_dist;
	public Genome gen;
	public int winSize;
	
	
	
	
	public GenomicRegionStats(Genome g, int win, String file_name, int min_d) throws IOException {
		this.gen = g;
		this.winSize = win;
		this.min_dist = min_d;
		BufferedReader reader = new BufferedReader(new FileReader(file_name));
	    String line;
	    while ((line = reader.readLine()) != null) {
	    	String[] pieces = line.split("\t");
	    	points_lists.add(RegionFileUtilities.loadPeaksFromPeakFile(gen, pieces[1], winSize));
	    	factor_string.add(pieces[0]);
	    }
	    reader.close();
	}
	
	public void setIntersectNumbers(String id_string){
		String[] idss = id_string.split(":");
		List<Integer> index_list = new ArrayList<Integer>();
		for(int i=0; i<idss.length; i++){
			index_list.add(Integer.parseInt(idss[i]));
		}
		String index_key= join(index_list,":");
		
		System.out.println(index_key+"\t"+this.getIntersectIndexes(index_list).size());
	}
	
	
	public String join(List<Integer> list, String conjunction)
	{
	   StringBuilder sb = new StringBuilder();
	   boolean first = true;
	   for (int item : list)
	   {
	      if (first)
	         first = false;
	      else
	         sb.append(conjunction);
	      sb.append(factor_string.get(item));
	   }
	   return sb.toString();
	}
	
	
	
	public List<Integer> getIntersectIndexes(List<Integer> list_ids){
		List<Integer> ret = new ArrayList<Integer>();
		List<Point> ref_list = this.points_lists.get(list_ids.get(0));
		for(int i=0; i<ref_list.size(); i++){
			boolean intersect = true;
			Point ref_point = ref_list.get(i);
			for(int l=1; l<list_ids.size(); l++){
				int mid_ind = this.getMinDistIndex(ref_point, this.points_lists.get(list_ids.get(l)));
				if(!ref_point.getChrom().equals(this.points_lists.get(list_ids.get(l)).get(mid_ind).getChrom())){
					intersect =false;
				}else{
					int dist = ref_point.distance(this.points_lists.get(list_ids.get(l)).get(mid_ind));
					if(dist >  this.min_dist){
						intersect = false;
					}
				}
			}
			if(intersect){
				ret.add(i);
			}
		}
		return ret;
	}
	
	public Integer getMinDistIndex(Point a, List<Point> pts){
		int index=0;
		int min = Integer.MAX_VALUE;
		for(int i=0; i<pts.size(); i++){
			if(a.getChrom().equals(pts.get(i).getChrom())){
				int dist = a.distance(pts.get(i));
				if(dist < min){
					min = dist;
					index = i;
				}	
			}
		}
		return index;
	}

	public static void main(String[] args) throws NotFoundException, IOException{
		ArgParser ap = new ArgParser(args);
		Genome g = null;
		if(ap.hasKey("species")){
			Pair<Species, Genome> pair = Args.parseGenome(args);
			if(pair != null){
				
				g = pair.cdr();
			}
		}else{
			if(ap.hasKey("geninfo") || ap.hasKey("g")){
				//Make fake genome... chr lengths provided
				String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
				g = new Genome("Genome", new File(fName), true);
			}else{
				g = null;
			}
		}
		
		int w= ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue() : 200;
		int mindd = ap.hasKey("mind") ? new Integer(ap.getKeyValue("mind")).intValue() : 100;
		String input_file = ap.getKeyValue("if");
		
		String intersect_ids = ap.getKeyValue("ids");
		
		GenomicRegionStats runner = new GenomicRegionStats(g,w,input_file, mindd);
		runner.setIntersectNumbers(intersect_ids);
		
		
		
	
	}
	
}
