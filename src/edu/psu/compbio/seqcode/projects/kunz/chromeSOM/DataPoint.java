package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;

import java.util.Vector;

public class DataPoint 
{
	public String look;
	public Node myNode;
	public MiniNode myMini;
	public boolean count;
	public int x,y;
	public double mag;
	public double[] g;
	public String name;
	public int chrome, minLocus, maxLocus;
	
	public DataPoint()
	{
		count = true;
	}
	public DataPoint(int chr, int min, int max, MiniNode n)
	{
		chrome = chr;
		minLocus = min;
		maxLocus = max;
		myMini = n;
	}
	public DataPoint(Node n)
	{
		myNode = n;
		count = true;
	}
	public DataPoint(String s)
	{
		name = s;
		look = "";
	}
	public void updateMag()
	{
		double magD = 0;
		for(int i = 0; i<g.length; i++)
		{
			magD += (g[i])*(g[i]);
		}
		magD = Math.sqrt(magD);
		mag = magD;
	}
	public void gInit(int length) 
	{
		g = new double[length];
	}
	public void setName(String s)
	{
		name = s;
	}
}