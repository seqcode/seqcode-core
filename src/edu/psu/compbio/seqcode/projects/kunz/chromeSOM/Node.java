package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;
import java.awt.Color;
import java.awt.Polygon;
import java.util.ArrayList;
import java.util.Vector;

public class Node
{
	public int xLoc, yLoc, xVal, yVal;
	public double mag;
	public ArrayList<ArrayList<Node>> neighborHoods;
	public ArrayList<Node> neighbors;
	public ArrayList<Node> neighbors2;
	public ArrayList<Node> neighbors3;
	public ArrayList<DataPoint> dataPoints;
	public ArrayList<DataPoint> countedPoints;
	public Polygon p;
	public Color color;
	public double[] g;
	public Node(int xa, int ya)
	{
		neighborHoods = new  ArrayList<ArrayList<Node>>();
		dataPoints = new ArrayList<DataPoint>();
		countedPoints = new ArrayList<DataPoint>();
		xVal = xa;
		yVal = ya;
		neighbors = new ArrayList<Node>();
		neighbors2 = new ArrayList<Node>();
		neighbors3 = new ArrayList<Node>();
		color = new Color((int)(Math.random()*255),(int)(Math.random()*255),(int)(Math.random()*255));
	}
	public void initialize(DataPoint d)
	{
		g = d.g;
	}
	public void polygonMaker(int[] x, int[] y, int i)
	{
		p = new Polygon(x,y,i);
	}
	public void setColor(Color c)
	{
		color = c;
	}
	public void setValues(Node n)
	{
		xLoc = n.xLoc;
		yLoc = n.yLoc;
		xVal = n.xVal;
		yVal = n.yVal;
		g = n.g;
	}
	public void neighborHoods() 
	{
		ArrayList<Node> me = new ArrayList<Node>();
		me.add(this);
		neighborHoods.add(me);
		neighborHoods.add(neighbors);
		neighborHoods.add(neighbors2);
		neighborHoods.add(neighbors3);
	}
	public boolean neighborHoodContains(Node add)
	{
		boolean bb = false;
		for(int i = 0; i < neighborHoods.size(); i++)
		{
			for(int j = 0; j< neighborHoods.get(i).size(); j++)
			{
				if(neighborHoods.get(i).indexOf(add)!=-1)
					bb = true;
			}
		}
		return bb;
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
}