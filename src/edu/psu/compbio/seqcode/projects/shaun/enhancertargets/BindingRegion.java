package edu.psu.compbio.seqcode.projects.shaun.enhancertargets;

import java.util.ArrayList;

import edu.psu.compbio.seqcode.genome.location.NamedPoint;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;

public class BindingRegion implements Comparable<BindingRegion>{

	public Region reg;
	public Point midpoint;
	public ArrayList<NamedPoint> points = new ArrayList<NamedPoint>();
	
	public BindingRegion(NamedPoint p){
		midpoint=p;
		reg = new Region(p.getGenome(), p.getChrom(), p.getLocation()-1, p.getLocation()+1);
		points.add(p);
	}
	
	//Add new point -- extend the region and update midpoint
	public void addPoint(NamedPoint p){
		if(reg.getChrom().equals(p.getChrom())){
			points.add(p);
			int min=p.getLocation(), max=p.getLocation();
			int sum = p.getLocation(), count = 1;
			for(NamedPoint x : points){
				if(x.getLocation()<p.getLocation())
					min=x.getLocation();
				else if(x.getLocation()>p.getLocation())
					max = x.getLocation();
				sum += x.getLocation();
				count++;
			}
			reg = new Region(reg.getGenome(), reg.getChrom(), min-1, max+1);
			midpoint = new Point(reg.getGenome(), reg.getChrom(), sum/count);
		}
	}
	//Add new points
	public void addPoint(BindingRegion br){
		for(NamedPoint x : br.points){
			addPoint(x);
		}
	}
	
	//String representation
	public String toString(){
		String res = reg.getLocationString()+"\t"+midpoint.getLocationString()+"\t"+points.size();
		for(NamedPoint np : points){
			res = res+"\t"+np.getName()+":"+np.getLocationString();
		}
		return res;
	}
	
	//Rank according to midpoints 
	public int compareTo(BindingRegion e) {
		return(midpoint.compareTo(e.midpoint));
	}
}
