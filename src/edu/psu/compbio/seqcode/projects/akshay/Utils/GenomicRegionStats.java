package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

public class GenomicRegionStats {
	public List<List<Point>> points_lists = new ArrayList<List<Point>>();
	public List<String> factor_string = new ArrayList<String>();
	public int min_dist;
	public Genome gen;
	public int winSize;
	
	public Map<String, Integer> interesect_numbers = new HashMap<String, Integer>();
	
	
	public GenomicRegionStats(Genome g, int win, String file_name, int min_d) throws IOException {
		this.gen = g;
		this.winSize = win;
		this.min_dist = min_d;
		BufferedReader reader = new BufferedReader(new FileReader(file_name));
	    String line;
	    while ((line = reader.readLine()) != null) {
	    	String[] pieces = line.split("\t");
	    	points_lists.add(Utilities.loadPeaksFromPeakFile(gen, pieces[1], winSize));
	    	factor_string.add(pieces[0]);
	    }
	    reader.close();
	}
	
	
	public void setIntersectNumbers(){
		for(int i=0; i< points_lists.size()-1; i++){
			for(int j=i; j<points_lists.size(); j++){
				List<Integer> index_list = this.getIndexString(i, j);
				if(i==j){
					index_list.add(i);
				}
				String index_key= join(index_list,":");
				this.interesect_numbers.put(index_key, this.getIntersectIndexes(index_list).size());
				
			}
		}
	}
	
	public void printOut(){
		for(String s : this.interesect_numbers.keySet()){
			System.out.println(s+"\t"+Integer.toString(this.interesect_numbers.get(s)));
		}
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
	
	private List<Integer> getIndexString(int start, int rem){
		List<Integer> ret = new ArrayList<Integer>();
		for(int i=start; i<this.points_lists.size(); i++){
			if(i != rem){
				ret.add(i);
			}
		}
		
		return ret;
	}
	
	public List<Integer> getIntersectIndexes(List<Integer> list_ids){
		List<Integer> ret = new ArrayList<Integer>();
		List<Point> ref_list = this.points_lists.get(list_ids.get(0));
		for(int i=0; i<ref_list.size(); i++){
			boolean intersect = true;
			Point ref_point = ref_list.get(i);
			for(int l=1; l<list_ids.size(); l++){
				int mid_ind = this.getMinDistIndex(ref_point, this.points_lists.get(l));
				if(!ref_point.getChrom().equals(this.points_lists.get(l).get(mid_ind).getChrom())){
					intersect =false;
				}else{
					int dist = ref_point.distance(this.points_lists.get(l).get(mid_ind));
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
			Pair<Organism, Genome> pair = Args.parseGenome(args);
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
		
		GenomicRegionStats runner = new GenomicRegionStats(g,w,input_file, mindd);
		runner.setIntersectNumbers();
		runner.printOut();
		
		
	
	}
	
}
