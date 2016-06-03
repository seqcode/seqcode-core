package org.seqcode.projects.kunz.isotop;


import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

public class Reader 
{
	public int count;
	public ArrayList<Prototype> protos;
	public String file;
	public String lookup;
	
	public Reader(String fileName)//, String look)
	{
		file = fileName;
		//lookup = look;
		protos = new ArrayList<Prototype>();
		count = 0;
	}
	public ArrayList<Prototype> matrixReader()
	{
		BufferedReader br = null;
		
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
            line = br.readLine();
	        while (line != null) 
	        {
	        	line = br.readLine();
	            makePoint(line);
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
	    return protos;
	}
	public void makePoint(String s)
	{
		if(s !=null && !s.substring(0, s.indexOf("\t")).isEmpty())
		{
			count++;
			//System.out.println(s);
			String name = s.substring(0,s.indexOf("\t"));
			s = s.substring(s.indexOf("\t")+1,s.length());
			
			String[] temp = s.split("\\t");
			double[] funk = new double[temp.length];
			try {
				for(int t = 0; t< temp.length; t++)
				{
					funk[t] = ((double)Double.parseDouble(temp[t]));
				}
				Prototype p = new Prototype(funk,name);
				protos.add(p);
				
			} catch (NumberFormatException e) {
			}	
		}
	}
}
