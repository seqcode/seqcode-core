package edu.psu.compbio.seqcode.projects.shaun;

import java.util.ArrayList;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;

public class ChromBlock{
	public ChromBlock(Region r){region=r;}
	public Region region;
	public int numPoints=0;
	public HashMap<Point, Integer> pointIndex = new HashMap<Point, Integer>();
	public ArrayList<Point> points=new ArrayList<Point>();
	public ArrayList<ArrayList<Double>> data = new ArrayList<ArrayList<Double>>();
}