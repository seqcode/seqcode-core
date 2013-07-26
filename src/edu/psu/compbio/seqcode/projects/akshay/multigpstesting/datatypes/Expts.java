package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

import java.io.*;
import java.util.*;


public class Expts {
	private Map<Point,String> isolatedPoints= new HashMap<Point,String>();
	private Map<Point,String> listOfPoints = new HashMap<Point,String>();
	private Map<Range,List<Point>> rangeMap = new HashMap<Range,List<Point>>();
	private List<Range> blacklist = new ArrayList<Range>();
	
	/* Fills rangeMap and list of points */
	public void mapPointstoRanges(String Filename) throws IOException{
		BufferedReader br = null;
		String currentline;
		br = new BufferedReader(new FileReader(Filename));
		currentline = br.readLine();
		while(currentline != null){
			if(!currentline.startsWith("#")){
				Point currentPoint = new Point(currentline.split("\t")[0]);
				listOfPoints.put(currentPoint, "");
				if(rangeMap.containsKey(currentPoint.overlaps.getFirst())){
					rangeMap.get(currentPoint.overlaps.getFirst()).add(currentPoint);
				}
				else{
					List<Point> temp = new LinkedList<Point>(Arrays.asList(currentPoint));
					rangeMap.put(currentPoint.overlaps.getFirst(), temp);
				}
				if(rangeMap.containsKey(currentPoint.overlaps.getLast())){
					rangeMap.get(currentPoint.overlaps.getLast()).add(currentPoint);
				}
				else{
					List<Point> temp = new LinkedList<Point>(Arrays.asList(currentPoint));
					rangeMap.put(currentPoint.overlaps.getLast(), temp);
				}
			}
			currentline = br.readLine();
		}
		br.close();
	}
	
	
	/* returns the list of isolated points */
	public Map<Point,String> getIsolatedPointMap(){
		return this.isolatedPoints;
	}
	
	public Map<Point,String> getListofPointMap(){
		return this.listOfPoints;
	}

	
	/* given a point, this function returns the distance to the point that is the closest in the list of isolated points*/
	public int getNearestDistance(Point currentPoint){
		int temp = Point.infdistance;
		for(Point iterPoint: isolatedPoints.keySet()){
			if(currentPoint.getDistancetoPoint(iterPoint) < temp){
				temp = currentPoint.getDistancetoPoint(iterPoint);
			}
		}
		return temp;
	}
	
	/* Blacklist here means points that are less than 500 bp apart form each other. This functions fills all the ranges that are blacklisted */
	public void fillBlacklist(){
		for (Range iterRange: rangeMap.keySet()){
			if(rangeMap.get(iterRange).size() > 1){
				blacklist.add(iterRange);
			}
		}
	}
	
	/* Based on the blacklist this function filters out all points that are not in the blacklist region and fills the isolatedpoints hashmap*/
	public void fillIsolatedPoints(Map<Range,Integer> unionBlacklist){
		this.isolatedPoints = this.listOfPoints;
		for(Point iterPoints: listOfPoints.keySet()){
			for(Range iterRange: unionBlacklist.keySet()){
				if(iterRange.includes(iterPoints)){
					isolatedPoints.remove(iterPoints);
					break;
				}
			}
		}
	}
	
	/* Returns the list of ranges that are blacklisted */
	public List<Range> getBlacklist(){
		return blacklist;
	}
	
	/* returns the isolatedpoints hashmap */
	public int getNumeberofIsolatedPoints(){
		return isolatedPoints.keySet().size();
	}
	
	/* returns the total number of points in this expt */
	public int getTotalNumberOfPoints(){
		return listOfPoints.size();
	}
	

}




class Tuple{
	public Range first;
	public Range last;
	Tuple(Range first, Range last){
		this.first = first;
		this.last = last;
	}
	
	
	Range getFirst(){
		return this.first;
	}
	Range getLast(){
		return this.last;
	}
	
	
}