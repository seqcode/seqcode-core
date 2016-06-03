package org.seqcode.projects.kunz.isotop;
import java.util.ArrayList;

public class Prototype 
{
	public String name;
	public String[] names;
	public int myIndex;
	public double pmag, weight, count;
	public double[] P_hiDVec;
	public double[] m;
	public int[] neighborIndexes;
	public ArrayList<Prototype> neighbors;
	public double[] neighborDists;
	double[][] loci;
	public Prototype(double[] vect, String n)
	{
		name = n;
		P_hiDVec = vect;
		m = new double[]{0,0};
		neighborIndexes = null;
		neighborDists = null;
		neighbors = new ArrayList<Prototype>();
		updateMag();
	}
	public Prototype(String s)
	{
		//System.out.println(s);
		ArrayList<String> naming = new ArrayList<String>();
		name = s.substring(0, s.indexOf("["));
		while(s.contains("-- "))
		{
			naming.add(s.substring(0,s.indexOf("-- ")));
			s = s.substring(s.indexOf("-- ")+3, s.length());
		}

		
		m = new double[2];
		m[0] = Double.parseDouble(s.substring(s.indexOf("[")+1,s.indexOf(";"))); 
		m[1] = Double.parseDouble(s.substring(s.indexOf("; ")+2,s.indexOf("]")));
		String[] namer = new String[naming.size()];
		for(int i = 0; i< naming.size(); i++)
		{
			namer[i] = naming.get(i);
		}
		setNames(namer);
				
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
		if(pmag == 0)
			pmag = .0000000001; //zero mag -> NaN for cosine similarity calculation
	}
	public void setName(String s)
	{
		name = s;
	}
	public void setNames(String[] s)
	{
		int co = 0;
		names = new String[s.length];
		for(int i =0; i< s.length; i++)
		{
			co++;
			names[i] = s[i];
		}
		loci = new double[co][3] ;
		for(int i =0; i < co; i++)
		{
			loci[i][0] = Integer.parseInt(s[i].substring(s[i].indexOf("chr")+3, s[i].indexOf(":")));;
			loci[i][1] = Integer.parseInt(s[i].substring(s[i].indexOf(":")+1, s[i].indexOf("-")));;
			loci[i][2] = Integer.parseInt(s[i].substring(s[i].indexOf("-")+1, s[i].length()));;
		}
		
	}
}
