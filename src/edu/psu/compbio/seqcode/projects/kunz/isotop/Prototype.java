package edu.psu.compbio.seqcode.projects.kunz.isotop;
import java.util.ArrayList;

public class Prototype 
{
	public String name;
	public String[] names;
	public int myIndex;
	public double pmag;
	public double[] P_hiDVec;
	public double[] m;
	public int[] neighborIndexes;
	public ArrayList<Prototype> neighbors;
	public double[] neighborDists;
	public Prototype(double[] vect, String n)
	{
		name =n;
		P_hiDVec = vect;
		m = new double[]{0,0};
		neighborIndexes = null;
		neighborDists = null;
		neighbors = new ArrayList<Prototype>();
		updateMag();
		if(pmag == 0)
			pmag = .00000001; //zero mag -> NaN for cosine similarity calculation
	}
	public Prototype(String s)
	{
		name = s.substring(0, s.indexOf("["));
		m = new double[2];
		m[0] = Double.parseDouble(s.substring(s.indexOf("[")+1,s.indexOf(";"))); 
		m[1] = Double.parseDouble(s.substring(s.indexOf("; ")+2,s.indexOf("]")));
				
	}
	public void updateMag()
	{
		double magD = 0;
		for(int i = 0; i<P_hiDVec.length; i++)
		{
			magD += (P_hiDVec[i])*(P_hiDVec[i]);
		}
		magD = Math.sqrt(magD);
		pmag = magD;
	}
	public void setName(String s)
	{
		name = s;
	}
}
