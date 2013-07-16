package edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes;

import java.util.*;


public class Expts {
	private Map<Point,String> isolatedPoints= new HashMap<Point,String>();
	private Map<Point,String> listOfPoints = new HashMap<Point,String>();
	private Map<Range,List<Point>> rangeMap = new HashMap<Range,List<Point>>();
	private List<Range> blacklist = new ArrayList<Range>();
	private String exptname;
	
	public void mapPointstoRnages(Point currentPoint, String exptname){
		listOfPoints.put(currentPoint, exptname);
		this.exptname = exptname;
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
	
	public String getExptname(){
		return this.exptname;
	}
	
	public int getNearestDistance(Point currentPoint){
		int temp = currentPoint.infdistance;
		for(Point iterPoint: isolatedPoints.keySet()){
			if(currentPoint.getDistancetoPoint(iterPoint) < temp){
				temp = currentPoint.getDistancetoPoint(iterPoint);
			}
		}
		return temp;
	}
	
	public void fillBlacklist(){
		for (Range iterRange: rangeMap.keySet()){
			if(rangeMap.get(iterRange).size() > 1){
				blacklist.add(iterRange);
			}
		}
	}
	
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
	
	public List<Range> getBlacklist(){
		return blacklist;
	}
	
	public int getNumeberofIsolatedPoints(){
		return isolatedPoints.keySet().size();
	}
	
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