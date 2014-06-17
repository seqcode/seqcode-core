package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;

public class GenomicRegionStats {
	public List<List<Point>> points_lists = new ArrayList<List<Point>>();
	public List<String> factor_string = new ArrayList<String>();
	public int min_dist;
	
	public Map<String, Integer> interesect_numbers = new HashMap<String, Integer>();
	
	public void setIntersectNumbers(){
		for(int i=0; i< points_lists.size(); i++){
			
		}
	}
	
	
	
	
	
	
	public List<Integer> getIntersectIndexes(List<Integer> list_ids){
		List<Integer> ret = new ArrayList<Integer>();
		List<Point> ref_list = this.points_lists.get(list_ids.get(0));
		for(int i=0; i<ref_list.size(); i++){
			boolean intersect = true;
			Point ref_point = ref_list.get(i);
			for(int l=1; l<list_ids.size(); l++){
				int mid_ind = this.getMinDistIndex(ref_point, this.points_lists.get(l));
				int dist = ref_point.distance(this.points_lists.get(l).get(mid_ind));
				if(dist >  this.min_dist){
					intersect = false;
				}
			}
			if(intersect){
				ret.add(i);
			}
		}
		return ret;
	}
	
	public Integer getMinDistIndex(Point a, List<Point> points_list){
		int index=0;
		int min = Integer.MAX_VALUE;
		for(int i=0; i<points_list.size(); i++){
			int dist = a.distance(points_list.get(i));
			if(dist < min){
				min = dist;
				index = i;
			}
		}
		return index;
	}

}
