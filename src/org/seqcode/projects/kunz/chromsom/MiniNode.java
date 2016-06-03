package org.seqcode.projects.kunz.chromsom;

import java.awt.Color;
import java.awt.Polygon;
import java.util.ArrayList;

public class MiniNode 
{
	public Polygon p;
	public ArrayList<DataPoint> dataPoints;
	public int xCoord, yCoord, xLoc, yLoc;
	public ArrayList<String> fullDPs;
	public ArrayList<DataPoint> counting;
	public ArrayList<DataPoint> notCounting;
	public Color color;
	public double weight;
	MiniNode(String s)
	{
		weight = 0;
		dataPoints = new ArrayList<DataPoint>();
		counting = new ArrayList<DataPoint>();
		notCounting = new ArrayList<DataPoint>();
		xCoord = Integer.parseInt(s.substring(s.indexOf("[")+1,s.indexOf(",")));
		yCoord = Integer.parseInt(s.substring(s.indexOf(",")+1,s.indexOf("]")));
		//System.out.println(xCoord +", "+yCoord);
	}
	public void addDP(String s)
	{
		if(s.contains(","))
			s = s.substring(0, s.indexOf(","));
		if (s.contains(")"))
			s = s.substring(0, s.indexOf(")"));
		if(s.contains("("))
			s = s.substring(s.indexOf("(")+1, s.length());
		if(s.contains(":")&&s.contains("-"))
		{
			int chrome = Integer.parseInt(s.substring(s.indexOf("chr")+3, s.indexOf(":")));
			int minLocus = Integer.parseInt(s.substring(s.indexOf(":")+1, s.indexOf("-")));
			int maxLocus = Integer.parseInt(s.substring(s.indexOf("-")+1, s.length()));
			{
				DataPoint d = new DataPoint(chrome,minLocus,maxLocus,this);
				d.name = s;
				dataPoints.add(d);
			}
		}
		
	}
	public void polygonMaker(int[] x, int[] y, int i)
	{
		p = new Polygon(x,y,i);
	}
}