package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;

import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

public class Matrixer 
{	
	int binsize, n;
	public int count;
	public ArrayList<Interaction> inter;
	public String file, file2, lander;
	public String lookup;
	int[] a, a1;
	int[][] bigpapa;
	
	public Matrixer()
	{
		binsize = 20000;
		String ffs = System.getProperty("user.dir")+"/MatrixLanding.txt";
		lander = ffs;
		String ff = System.getProperty("user.dir")+"/YeastInteractions-intra.txt";
		file = ff;
		String f = System.getProperty("user.dir")+"/YeastInteractions-inter.txt";
		file2 = f;
		inter = new ArrayList<Interaction>();
		count = 0;
	}
	public void matrixReader()
	{
		reader(file);
		reader(file2);
	    binGenome();
	}
	public void reader(String f)
	{
		BufferedReader br = null;
		
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(f));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
	        line = br.readLine();
	        while (line != null) 
	        {
	        	StringBuilder sb = new StringBuilder();
	            sb.append(line);
	            sb.append(System.lineSeparator());
	            line = br.readLine();
	            makePoint(sb.toString());
	        }
	    } catch (IOException e)
	    {
			e.printStackTrace();
		} finally 
		{if(br!=null)
			{
	        	try {
	        		br.close();
	        	} catch (IOException e) 
				{
				e.printStackTrace();
				}
			}
		}
	}
	public void binGenome()
	{
		n = findChromeNum();
		a = maxForEach(binsize,n);
		a1 = new int[n];
		System.out.println(n);
		int totalBins = 0;
		for(int i =0; i<a.length;i++)
		{
			a1[i] = totalBins;
			totalBins+=a[i];
		}
		for(int i =0; i<a.length;i++)
			System.out.println(a[i]);
		System.out.println(totalBins);
		bigpapa = new int[totalBins][totalBins];
		for(int i =0; i< inter.size();i++)
		{
			int[] b = findBin(inter.get(i));
			bigpapa[b[0]][b[1]] += (int)inter.get(i).g[4];
			bigpapa[b[1]][b[0]] += (int)inter.get(i).g[4];
		}
		writeFile(totalBins);
	}
	public int[] findBin(Interaction in)
	{
		int[] i ={0,0};
		
		i[0] = a1[(int)in.g[0]-1]+((int)in.g[1]/binsize);
		i[1] = a1[(int)in.g[2]-1]+((int)in.g[3]/binsize);
		return i;
	}
	public int findChromeNum()
	{
		int m=0;
		for(int i =0; i< inter.size();i++)
		{
			if(inter.get(i).g[0]>m)
				m = (int) inter.get(i).g[0];
			if(inter.get(i).g[2]>m)
				m = (int) inter.get(i).g[2];
		}
		return m;
	}
	public int[] maxForEach(int binsize, int m)
	{
		int[] r = new int[m];
		for(int i = 1; i<=m; i++)
		{
			int big = 0;
			for(int j = 0; j< inter.size(); j++)
			{
				if(inter.get(j).g[2] == i && inter.get(j).g[3]>big)
					big = (int)inter.get(j).g[3];
			}
			r[i-1]=big;
		}
		for(int i =0; i<r.length;i++)
		{
			r[i]= r[i]+(binsize-(r[i]%binsize));
			r[i]= r[i]/binsize;
		}
		return r;
	}
	public void makePoint(String s)
	{
		//System.out.println(s);
		count++;
		Interaction p = new Interaction();
		
		String[] temp = s.split("\\t");
		try {
			for(int t = 0; t< temp.length; t++)
			{
				p.g[t] = ((double)Double.parseDouble(temp[t]));
			}
			inter.add(p);
		} catch (NumberFormatException e) {
		}
	}
	public void writeFile(int totalBins)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff = new FileWriter(lander,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			printer.print("\t");
			for(int i = 0; i<n; i++)
			{
				for(int j = 0; j<a[i]; j++)
				{
					printer.print("chr" + (i+1) +":"+ j*binsize + "-" + (((j+1)*binsize)-1)  + "\t");
				}
			}
			printer.print("\n");
			int p =0 ; int o = 0;
			for(int i = 0; i<totalBins; i++)
			{
				printer.print("chr" + (p+1) +":"+ o*binsize + "-" + (((o+1)*binsize)-1)  + "\t");
				for(int j = 0; j<totalBins; j++)
				{
					printer.print(bigpapa[i][j] + "\t");
				}
				 o++;
				if(a[p]<=o)
				{o = 0; p++;}
				printer.print("\n");
			}
			printer.close();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
			System.out.print("no way");
		}
	}
	public static void main(String[] args)
	{
		Matrixer m = new Matrixer();
		m.matrixReader();
		/*for(int i =0; i< m.inter.size();i++)
		{
			for(int j = 0; j< m.inter.get(i).g.length;j++)
			{
				System.out.print(m.inter.get(i).g[j]+"\t");
			}
			System.out.print("\n");
		}*/
	}
	private class Interaction
	{
		private double[] g;
		private Interaction()
		{
			g = new double[6];
		}
	}
}
